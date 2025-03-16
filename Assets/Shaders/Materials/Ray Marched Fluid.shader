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
            
            // Allowed floating point inaccuracy
            #define EPSILON 1e-5
            
            struct Attributes
            {
                float4 positionOS : POSITION;
            };
            
            struct Varyings
            {
                float4 positionCS : SV_POSITION;
                float3 positionWS: TEXCOORD0;
            };

            struct Ray
            {
                float3 origin;
                float3 direction;
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

            // Constant
            const float AIR_IOR = 1.0003f;

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

            // Given a position in OS, returns the closest face normal
            float3 closest_face_normal(float3 positionOS)
            {
                // Computing the absolute value of each axis to find the most extended one
                float absolute_x = abs(positionOS.x);
                float absolute_y = abs(positionOS.y);
                float absolute_z = abs(positionOS.z);
                
                // Case in which X axis of the position is the largest
                if (absolute_x >= absolute_y && absolute_x >= absolute_z)
                {
                    // Returning the normal of the closest face perpendicular to the X axis
                    return float3(sign(positionOS.x), 0, 0);
                }
                
                // Case in which Y axis of the position is the largest
                if (absolute_y >= absolute_x && absolute_y >= absolute_z)
                {
                    // Returning the normal of the closest face perpendicular to the Y axis
                    return float3(0, sign(positionOS.y), 0);
                }

                // Returning the normal of the closest face perpendicular to the Z axis
                return float3(0, 0, sign(positionOS.z));
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
                float3 face_normal = closest_face_normal(positionOS);
                
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
            float total_density_along_rayOS (Ray ray, float step_size)
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
            
            // Fragment shader
            float4 frag (Varyings input) : SV_Target
            {
                // Computing the step sizes
                float step_size = sqrt(3.0f) / steps_amount;

                // Intializing the view ray
                Ray view_ray;
                view_ray.origin = TransformWorldToObject(input.positionWS);
                view_ray.direction = TransformWorldToObjectDir(normalize(input.positionWS - GetCameraPositionWS()));

                // Initializing the current position
                float3 current_positionOS = view_ray.origin + 1e-4 * view_ray.direction;

                // Initializing the color
                float4 final_color = 0;
                bool is_submerged = false;
                int current_surface_collisions = 0;
                
                // Raymarching to verify if the view ray intersect the fluid
                [loop]
                for(int i = 0; i < steps_amount; i++)
                {
                    // If I reached the boundaries of the cube, break
                    if(any(abs(current_positionOS) >= 0.5f + EPSILON))
                        break;
                    
                    // Sampling the density at the current position
                    float sampled_density = sample_density(current_positionOS);

                    // Verify if I reached a new surface
                    if(sampled_density >= 5e-1 && !is_submerged || sampled_density < 5e-1 && is_submerged)
                    {
                        // Increasing amount of surface collisions
                        current_surface_collisions++;
                        
                        // Computing the normal at the current point
                        float3 surface_normal = compute_surface_normal(current_positionOS);
                        
                        // Verifying if the direction and normal are correctly aligned
                        float3 correct_normal = dot(view_ray.direction, surface_normal) < 0 ?
                            surface_normal : -surface_normal;

                        // Verifying the index of refraction between the mediums
                        float correct_IOR = dot(view_ray.direction, surface_normal) < 0 ?
                            AIR_IOR / refraction_index : refraction_index / AIR_IOR;

                        // Computing the refraction direction using the current direction and the surface normal
                        float3 refracted_direction = refract(view_ray.direction, correct_normal, correct_IOR);

                        // Computing the reflection direction using the current direction and the surface normal
                        float3 reflected_direction = reflect(view_ray.direction, correct_normal);

                        // If the refraction direction is (0, 0, 0), refraction is impossible
                        if(all(refracted_direction == 0))
                        {
                            // Changing the current direction to the reflected direction
                            view_ray.direction = reflected_direction;
                            // Changing the final color to blue
                            final_color = float4(0, 0, 0.75, 1);
                        }
                        else
                        {
                            // Changing the current direction to the refracted direction
                            view_ray.direction = refracted_direction;
                            // Changing the final color to green
                            final_color = float4(0, 0.75, 0, 1);
                        }
                        
                        // If I reached the maximum surface collisions, break
                        if(current_surface_collisions >= max_surface_collisions)
                        {
                            break;
                        }
                        
                        // Flipping flag
                        is_submerged = !is_submerged;
                    }
                    
                    // Advancing the current position
                    current_positionOS += step_size * view_ray.direction;
                }

                return final_color;
            }
            
        ENDHLSL
        }
    }
}