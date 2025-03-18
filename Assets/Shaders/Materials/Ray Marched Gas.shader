Shader "Custom/Ray Marched Gas"
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
                    float sampled_density = sample_density(current_positionOS).x
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
