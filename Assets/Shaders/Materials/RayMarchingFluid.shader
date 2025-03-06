Shader "Custom/Ray Marching Fluid"
{
    Properties
    {
        // Density map input, provided by the fluid simulation
        _DensityMap ("Density Map", 3D) = "" {}
        // The ray marching step size
        _StepsAmount ("Steps Amount", Int) = 128
        // The cumulative alpha of the marching algorithm
        _DensityMultiplier ("Marching Alpha", Float) = 1.0
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
    
        Pass
        {
            Blend One OneMinusSrcAlpha
            
            HLSLPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma target 5.0

            // #include "UnityCG.cginc"
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            
            // Maximum number of raymarching samples
            #define MAX_STEP_COUNT 300
            #define SUN_STEPS_COUNT 75

            // Allowed floating point inaccuracy
            #define EPSILON 1e-5

            struct Attributes
            {
                float4 positionOS : POSITION;
                float3 normal : NORMAL;
                float2 uv : TEXCOORD0;
            };
            
            struct Varyings
            {
                float4 positionCS : SV_POSITION;
                float3 normal : NORMAL;
                float2 uv : TEXCOORD0;
                float3 positionOS: TEXCOORD1;
                float3 ray_directionOS : TEXCOORD2;
            };

            // Properties
            TEXTURE3D(_DensityMap);
            SAMPLER(sampler_DensityMap);
            int _StepsAmount;
            float _DensityMultiplier;

            // Variables
            float step_size;
            float3 scattering_coefficients;

            // Auxiliary functions
            float calculate_density_along_ray (float3 position, float3 direction, float step_size)
            {
                // Initializing the auxiliary variables
                float total_density = 0;

                // Raymarching through object space
                [loop]
                for(int i = 0; i < _StepsAmount; i++)
                {
                    // Computing the current position
                    float3 current_position = position + step_size * i * direction;

                    // If the current position left the cube boundaries, break
                    if(max(abs(current_position.x), max(abs(current_position.y), abs(current_position.z)))
                        > 0.5f + EPSILON)
                        break;

                    // Sampling the texture at the current position mapped to UV range [0,1]
                    float sampled_density = SAMPLE_TEXTURE3D(_DensityMap, sampler_DensityMap, current_position + 0.5).x
                        * step_size * _DensityMultiplier;

                    // Adding the sampled density to the total density
                    total_density += sampled_density;
                }

                return total_density;
            }
            
            // Vertex Shader
            Varyings vert (Attributes input)
            {
                // Initializing the ouput
                Varyings output;
                
                // Computing the vertex position in clip space
                output.positionCS = TransformObjectToHClip(input.positionOS.xyz);
                
                // Storing the vertex position in object space
                output.positionOS = input.positionOS;
                
                // Transforming the vetex in world space
                float3 positionWS = mul(UNITY_MATRIX_M, input.positionOS).xyz;

                // Computing the ray in world space
                float3 ray_directionWS = normalize(positionWS - _WorldSpaceCameraPos);
                
                // Computing the raay direction in world space
                output.ray_directionOS = normalize(mul(UNITY_MATRIX_I_M, float4(ray_directionWS, 0.0)));
                
                // Storing the uv coordinates
                output.uv = input.uv;
                
                // Storing the normal
                output.normal = input.normal;
                
                return output;
            }

            // Fragment shader
            float4 frag(Varyings input) : SV_Target
            {
                // Initializing the auxiliary variables
                float total_density = 0;
                float3 total_radiance = 0;
                
                // Raymarching through object space
                [loop]
                for(int i = 0; i < _StepsAmount; i++)
                {
                    // Computing the current position
                    float3 current_position = input.positionOS + step_size * i * input.ray_directionOS;
                    
                    // If the current position left the cube boundaries, break
                    if(max(abs(current_position.x), max(abs(current_position.y), abs(current_position.z)))
                        > 0.5f + EPSILON)
                        break;
                
                    // Sampling the texture at the current position mapped to UV range [0,1]
                    float sampled_density = SAMPLE_TEXTURE3D(_DensityMap, sampler_DensityMap, current_position + 0.5).x
                        * step_size * _DensityMultiplier;
                
                    // Adding the sampled density to the total density
                    total_density += sampled_density;

                    // Computing a fixed sun direction
                    float3 sun_direction = normalize(float3(1, 1, 1));
                    
                    // Calculating the density along the ray reaching the sun
                    float density_along_sun_ray = calculate_density_along_ray(current_position, sun_direction,
                        1 / 25.0f);
                    
                    // Calculating the radiance emitted from the sun reaching the current point
                    float3 sun_radiance = exp(-density_along_sun_ray * scattering_coefficients);
                
                    // Computed an approximation of the light scattered towards the camera
                    float3 scattered_light = sun_radiance * sampled_density * scattering_coefficients;
                
                    // Computing the radiance reaching the camera using the accumulated density
                    float3 current_radiance = exp(-total_density * scattering_coefficients);
                    
                    // Adding the radiance reaching the camera from current position to the total radiance
                    total_radiance += scattered_light * current_radiance;
                }
                
                return float4(total_radiance, 1.0f);
            }
        
        ENDHLSL
        }
    }
}
