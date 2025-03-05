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
            
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma target 5.0

            #include "UnityCG.cginc"
            
            // Maximum number of raymarching samples
            #define MAX_STEP_COUNT 500

            // Allowed floating point inaccuracy
            #define EPSILON 1e-5
            
            struct v2f
            {
                float4 position : SV_POSITION;
                float3 normal : NORMAL;
                float2 uv : TEXCOORD0;
                float3 positionOS: TEXCOORD1;
                float3 ray_directionOS : TEXCOORD2;
            };

            // Proprties
            sampler3D _DensityMap;
            int _StepsAmount;
            float _DensityMultiplier;

            // Variables
            float step_size;
            
            // Vertex Shader
            v2f vert (appdata_base input)
            {
                // Initializing the ouput
                v2f output;
                
                // Computing the vertex position in clip space
                output.position = UnityObjectToClipPos(input.vertex);
                
                // Storing the vertex position in object space
                output.positionOS = input.vertex;
                
                // Transforming the vetex in world space
                float3 positionWS = mul(unity_ObjectToWorld, input.vertex).xyz;

                // Computing the ray in world space
                float3 ray_directionWS = normalize(positionWS - _WorldSpaceCameraPos);
                
                // Computing the raay direction in world space
                output.ray_directionOS = normalize(mul(unity_WorldToObject, float4(ray_directionWS, 0.0)));
                
                // Storing the uv coordinates
                output.uv = input.texcoord;
                
                // Storing the normal
                output.normal = input.normal;
                
                return output;
            }

            // Fragment shader
            float4 frag(v2f input): SV_Target
            {
                // Initializing the color
                float total_density = 0;
                float3 total_radiance = 0;

                // Initializing the distance travelled within the simulation boundaries
                

                // Raymarching through object space
                for(int i = 0; i < min(_StepsAmount, MAX_STEP_COUNT); i++)
                {
                    // Computing the current position
                    float3 current_position = input.positionOS + step_size * i * input.ray_directionOS;
                    
                    // If the current position left the cube boundaries, break
                    if(max(abs(current_position.x), max(abs(current_position.y), abs(current_position.z)))
                        > 0.5f + EPSILON)
                        break;

                    // Sampling the texture at the current position mapped to UV range [0,1]
                    float sampled_density = tex3D(_DensityMap, current_position + 0.5).x
                        * step_size * _DensityMultiplier;

                    // Adding the sampled density to the total density
                    total_density += sampled_density;

                    // Computed an approximation of the scattered light towards the camera
                    float3 scattered_light = float3(1, 1, 1) * sampled_density;

                    // Computing the radiance reaching the camera using the accumulated density
                    float3 current_radiance = exp(-total_density);

                    total_radiance += scattered_light * current_radiance;
                }
 
                return float4(total_radiance, 1.0f);
            }
        
        ENDCG
        }

    }
}
