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
            int internal_steps_amount;
            float density_multiplier;
            float3 scattering_coefficients;
            float3 sun_direction;

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
                return SAMPLE_TEXTURE3D(_DensityMap, sampler_DensityMap, positionOS + 0.5f);
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
            float3 surface_normal(float3 positionOS)
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
                    float sampled_density = sample_density(current_position).x
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
                float internal_step_size = sqrt(3.0f) / internal_steps_amount;

                // Intializing the view ray
                Ray view_rayOS;
                view_rayOS.origin = TransformWorldToObject(input.positionWS);
                view_rayOS.direction = TransformWorldToObjectDir(normalize(input.positionWS - GetCameraPositionWS()));

                // Initializing the color
                float4 final_color = 0;
                
                // Raymarching to verify if the view ray intersect the fluid6
                [loop]
                for(int i = 0; i < steps_amount; i++)
                {
                    // Computing the current position
                    float3 current_positionOS = view_rayOS.origin + i * step_size * view_rayOS.direction;

                    // If I reached the boundaries of the cube, break
                    if(any(abs(current_positionOS) >= 0.5f + EPSILON))
                        break;
                    

                    // Sampling the density at the current position
                    float sampled_density = sample_density(current_positionOS).r;

                    // If the density is above the threshold, compute the normal of the surface
                    if(sampled_density >= 5e-1)
                    { 
                        final_color = float4((surface_normal(current_positionOS) + 1) * 0.5, 1);
                        break;
                    }
                }

                return final_color;
            }

            // Fragment shader
            float4 old_frag (Varyings input) : SV_Target
            {
                // Initializing the auxiliary variables
                float total_density = 0;
                float3 total_radiance = 0;

                // Computing the step sizes
                float step_size = sqrt(3.0f) / steps_amount;
                float internal_step_size = sqrt(3.0f) / internal_steps_amount;

                // Intializing the view ray
                Ray view_rayOS;
                view_rayOS.origin = TransformWorldToObject(input.positionWS);
                view_rayOS.direction = TransformWorldToObjectDir(normalize(input.positionWS - GetCameraPositionWS()));

                // Raymarching through object space
                [loop]
                for(int i = 0; i < steps_amount; i++)
                {
                    // Computing the current position
                    float3 current_positionOS = view_rayOS.origin + i * step_size * view_rayOS.direction;
                    
                    // If I reached the boundaries of the cube, break
                    if(any(abs(current_positionOS) >= 0.5f + EPSILON))
                        break;

                    // Sampling the density at the current position
                    float sampled_density = SAMPLE_TEXTURE3D(_DensityMap, sampler_DensityMap, current_positionOS + 0.5f).r
                        * density_multiplier * step_size;

                    // Adding the sampled density to the total
                    total_density += sampled_density;

                    // Initializing the sun ray reaching the current position in the fluid
                    Ray sun_rayOS;
                    sun_rayOS.origin = current_positionOS;
                    sun_rayOS.direction = sun_direction;

                    // Computing the density along the sun ray
                    float3 density_along_sun_ray = total_density_along_rayOS(sun_rayOS, internal_step_size);

                    // Computing an approximation of the radiance transmitted from the sun to the current point
                    float3 sun_radiance = exp(-density_along_sun_ray * scattering_coefficients);
                    
                    // Computing an approximation of the light scattered towards the camera
                    float3 scattered_light_towards_camera = sun_radiance * sampled_density * scattering_coefficients;

                    // Computing an approximation of the radiance transmitted towards the camera
                    float3 current_radiance = exp(- total_density * scattering_coefficients);

                    // Adding the radiance that reaches the camera from the current position
                    total_radiance += scattered_light_towards_camera * current_radiance;
                }

                return float4(total_radiance, 1.0f);
            }
        
        ENDHLSL
        }
    }
}
