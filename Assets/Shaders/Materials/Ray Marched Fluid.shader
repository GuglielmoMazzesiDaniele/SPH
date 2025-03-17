Shader "Custom/Ray Marched Fluid"
{
    Properties
    {
        // Density map input, provided by the fluid simulation
        _DensityMap ("Density Map", 3D) = "" {}
    }
    SubShader
    {
        Tags {"Queue"="Transparent" "RenderType"="Transparent"}
    
        Pass
        {
            Blend SrcAlpha OneMinusSrcAlpha
            
            HLSLPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma target 5.0

            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            #include "../Auxiliary/Auxiliary.hlsl"
            
            // Allowed floating point inaccuracy
            #define EPSILON 1e-5
            #define AIR_IOR 1.0003

            // Structs

            struct surface_collision
            {
                Ray reflected_ray;
                Ray refracted_ray;
                float reflected_intensity;
                float refracted_intensity;
                bool is_refraction_impossible;
            };
            
            struct Attributes
            {
                float4 positionOS : POSITION;
            };
            
            struct Varyings
            {
                float4 positionCS : SV_POSITION;
                float3 positionWS: TEXCOORD0;
            };

            // Properties
            TEXTURE3D(_DensityMap);
            SAMPLER(sampler_DensityMap);

            // Variables
            int steps_amount;
            int max_surface_collisions;
            float density_multiplier;
            float refraction_index;
            float3 scattering_coefficients;
            float4x4 floor_world_to_object;
            
            // Vertex Shader
            Varyings vert (Attributes input)
            {
                // Initializing the ouput
                Varyings output;
                
                // Computing the vertex position in clip space
                output.positionCS = TransformObjectToHClip(input.positionOS.xyz);

                // Storing the position in object space
                output.positionWS = TransformObjectToWorld(input.positionOS);
                
                return output;
            }
            
            // Given a position in OS, returns the density at the given position
            float sample_density(float3 positionOS)
            {
                return SAMPLE_TEXTURE3D(_DensityMap, sampler_DensityMap, positionOS + 0.5f).x;
            }

            // Given a ray in OS, sample light received from the enviroment
            float3 sample_enviroment(Ray ray)
            {
                // Converting the Ray in WS
                ray.origin = TransformObjectToWorld(ray.origin);
                ray.direction = normalize(TransformObjectToWorldDir(ray.direction));

                // Trying to sample the floor
                float2 floor_uv = ray_floor_intersection(ray, floor_world_to_object);

                // If the ray intersect the floor, returns it color
                if(all(floor_uv != -1))
                    return floor_uv_to_color(floor_uv);
                
                // TODO: Sampling the sky

                // Default color: black
                return float3(0, 0, 0);
            }

            // Given a position in OS, return its surface normal
            float3 compute_surface_normal(float3 positionOS)
            {
                // Auxiliary constants
                const float sampling_offset = 1e-3;
                const float max_boundaries_smoothing_distance = 8e-3;

                // Computing the offsets per axis
                float3 offsetX = float3(1, 0, 0) * sampling_offset;
                float3 offsetY = float3(0, 1, 0) * sampling_offset;
                float3 offsetZ = float3(0, 0, 1) * sampling_offset;

                // Computing the gradient along all axis
                float dx = sample_density(positionOS - offsetX)
                    - sample_density(positionOS + offsetX);
                float dy = sample_density(positionOS - offsetY)
                    - sample_density(positionOS + offsetY);
                float dz = sample_density(positionOS - offsetZ)
                    - sample_density(positionOS + offsetZ);

                // Normalizing the gradients to obtain the volume normal
                float3 volume_normal = normalize(float3(dx, dy, dz));

                // Obtaining the closest face normal to the current position
                float3 face_normal = cube_closest_face_normal(positionOS);
                
                // Computing the distance to boundaries
                float3 boundaries_distance = 0.5f - abs(positionOS);

                // Extracting the distance to the closest boundary
                float min_boundaries_distance = min(boundaries_distance.x,
                    min(boundaries_distance.y, boundaries_distance.z));

                // Computing the weight of the face normal based on the distance to it
                float face_weight = (1 - smoothstep(0, max_boundaries_smoothing_distance, min_boundaries_distance));

                // Computing the interpolated normal between volume and face normal
                float3 interpolated_normal = normalize(volume_normal * (1 - face_weight)
                    + face_normal * face_weight);
                
                return interpolated_normal;
            }

            // Given a ray in object space, computes the total fluid density along the ray
            float total_density_along_ray (Ray ray, float step_size)
            {
                // Initializing the auxiliary variables
                float total_density = 0;

                // Raymarching through object space
                [loop]
                for(int i = 0; i < steps_amount; i++)
                {
                    // Computing the current position
                    float3 current_position = ray.origin + step_size * i * ray.direction;

                    // If the current position left the cube boundaries, break
                    if(any(abs(current_position) > 0.5f + EPSILON))
                        break;

                    // Sampling the texture at the current position mapped to UV range [0,1]
                    float sampled_density = sample_density(current_position)
                        * step_size * density_multiplier;

                    // Adding the sampled density to the total density
                    total_density += sampled_density;
                }

                return total_density;
            }

            // Given a density, returns the transmitted light along the density
            float3 transmittance_decay (float density)
            {
                return exp(-density * scattering_coefficients);
            }

            // Given a normal and an incident direction in the same hemisphere, returns the fresnel's refraction
            // coefficient
            float schlick_approximation (float3 incident_direction, float3 normal)
            {
                // Intializing the F0 value, set to 0.02 for air <-> water interactions (equation is symmetric)
                const float f_0 = 0.02;
                
                // Computing the cosine of the angle between the vectors
                float cos_theta = abs(dot(normal, incident_direction));

                // Computing the fresnel effect coefficient
                return f_0 + (1.0 - f_0) * pow(1.0 - cos_theta, 5);
            }

            // Given a normal, an incident direction and a collision position, 
            // returns the struct containing the data regarding the surface collision
            surface_collision calculate_surface_relect_and_refract (float3 normal, float3 incident_direction,
                float3 collision_position)
            {
                // Verifying if the direction and normal are correctly aligned
                float3 correct_normal = dot(incident_direction, normal) < 0 ?
                    normal : -normal;

                // Verifying the index of refraction between the mediums
                float correct_IOR = dot(incident_direction, normal) < 0 ?
                    AIR_IOR / refraction_index : refraction_index / AIR_IOR;

                // Initializing the reflected ray
                Ray reflected_ray;
                reflected_ray.origin = collision_position;
                reflected_ray.direction = reflect(incident_direction, correct_normal);

                // Initializing the refracted ray
                Ray refracted_ray;
                refracted_ray.origin = collision_position;
                refracted_ray.direction = refract(incident_direction, correct_normal, correct_IOR);
                
                // Computing the reflected intensity using Schlick approximation
                float reflected_intensity = schlick_approximation(incident_direction, correct_normal);

                // Computing the refracted intensity
                float refracted_intensity = 1.0 - reflected_intensity;
                
                // Initializing the returned struct
                surface_collision collision_info;

                // Storing the computed data
                collision_info.reflected_ray = reflected_ray;
                collision_info.refracted_ray = refracted_ray;
                collision_info.reflected_intensity = reflected_intensity;
                collision_info.refracted_intensity = refracted_intensity;
                collision_info.is_refraction_impossible = all(refracted_ray.direction == 0);

                // Returning the computed data
                return collision_info;
            }
            
            // Fragment shader
            float4 frag (Varyings input) : SV_Target
            {
                // Computing the step sizes
                float step_size = sqrt(3.0f) / steps_amount;

                // Intializing the travelling direction
                float3 ray_marching_direction = TransformWorldToObjectDir(
                    normalize(input.positionWS - GetCameraPositionWS()));

                // Initializing the current position
                float3 ray_marching_position = TransformWorldToObject(input.positionWS) + 1e-4 * ray_marching_direction;

                // Auxiliary variables
                float3 final_color = 0;
                float collected_density = 0;
                float remaining_transmittance = 1;
                bool is_submerged = false;
                int current_surface_collisions = 0;
                
                // Raymarching within the simulation
                [loop]
                for(int i = 0; i < steps_amount; i++)
                {
                    // Sampling the density at the current position
                    float sampled_density = sample_density(ray_marching_position);

                    // Adding the sampled density to the accumulated density along current path
                    collected_density += sampled_density * density_multiplier * step_size;
     
                    // Verify if I reached the fluid surface, either from inside or outside of it
                    if((sampled_density >= 2.5e-1 && !is_submerged) || (sampled_density <= 2.5e-1 && is_submerged))
                    {
                        // Increasing amount of surface collisions
                        current_surface_collisions++;
                        
                        // Flipping flag
                        is_submerged = !is_submerged;

                        // Reducing the transmittance by the amount of density the ray passed through
                        remaining_transmittance *= transmittance_decay(collected_density);

                        // Resetting the collected density
                        collected_density = 0;
                         
                        // Computing the normal at the current point
                        float3 surface_normal = compute_surface_normal(ray_marching_position);

                        // If I reached the maximum surface collisions, break
                        if(current_surface_collisions >= max_surface_collisions)
                            break;

                        // Computing the surface collision info
                        surface_collision collision_info = calculate_surface_relect_and_refract(surface_normal,
                            ray_marching_direction, ray_marching_position);

                        // If refraction is not possible, follow the reflected path
                        if(collision_info.is_refraction_impossible)
                        {
                            // Changing the current direction to the reflected direction
                            ray_marching_direction = collision_info.reflected_ray.direction;
                        }
                        // If refraction is possible, approximated lesser path and follow greater one
                        else
                        {
                            // Computing the density along the reflected ray
                            float density_along_reflected = total_density_along_ray(collision_info.reflected_ray,
                                step_size * 10);
                            
                            // Computing the density along the
                            float density_along_refracted = total_density_along_ray(collision_info.refracted_ray,
                                step_size * 10);

                            // Deciding which path to follow based on the intensity of the path and the total fluid
                            // density it will traverse
                            bool is_refracted_path_better = density_along_refracted * collision_info.refracted_intensity
                                > density_along_reflected * collision_info.reflected_intensity;
                            
                            // If the refracted path is better, change the direction of the travelling ray to follow it.
                            // Otherwise change the direction to match the reflected path.
                            if(is_refracted_path_better)
                            {
                                // Changing direction to refracted
                                ray_marching_direction = collision_info.refracted_ray.direction;

                                // Approximating the contribution of the discated reflected path
                                final_color += sample_enviroment(collision_info.reflected_ray)
                                    * remaining_transmittance * transmittance_decay(density_along_reflected)
                                    * collision_info.reflected_intensity;

                                // Reducing the remaining trasmittance
                                remaining_transmittance *= collision_info.refracted_intensity;
                            }
                            else
                            {
                                // Changing direction to refracted
                                ray_marching_direction = collision_info.reflected_ray.direction;

                                // Approximating the contribution of the discated refracted path
                                final_color += sample_enviroment(collision_info.refracted_ray)
                                    * remaining_transmittance * transmittance_decay(density_along_reflected)
                                    * collision_info.refracted_intensity;

                                // Reducing the remaining transmittance
                                remaining_transmittance *= collision_info.reflected_intensity;
                            }
                        }
                    }
                    
                    // Advancing the current position
                    ray_marching_position += step_size * ray_marching_direction;
                    
                    // If I reached the boundaries of the cube, break
                    if(any(abs(ray_marching_position) >= 0.5f))
                        break;
                }

                return float4(final_color, 1);
            }
            
        ENDHLSL
        }
    }
}