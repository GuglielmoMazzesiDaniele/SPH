Shader "Custom/Ray Marching Fluid"
{
    Properties
    {
        // Density map input, provided by the fluid simulation
        _DensityMap ("Density Map", 3D) = "" {}
        // The ray marching step size
        _StepSize ("Step Size", Float) = 0.01
        // The cumulative alpha of the marching algorithm
        _DensityMultiplier ("Marching Alpha", Float) = 1.0
        // The minimum density value to depict
        _MinDensityValue ("Minumum Density Value", Float) = 150.0
        // The maximum density value to depict
        _MaxDensityValue ("Maximum Density Value", Float) = 1000.0
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
            #define MAX_STEP_COUNT 128

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

            // Texture and Uniforms
            sampler3D _DensityMap;
            float _StepSize;
            float _DensityMultiplier;
            float _MinDensityValue;
            float _MaxDensityValue;

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

            float4 BlendUnder(float4 color, float4 newColor)
            {
                color.rgb += (1.0 - color.a) * newColor.a * newColor.rgb;
                color.a += (1.0 - color.a) * newColor.a;
                return color;
            }

            float4 frag(v2f input): SV_Target
            {
                // Initializing the color
                float3 total_density = 0;

                // Raymarching through object space
                for(int i = 0; i < MAX_STEP_COUNT; i++)
                {
                    // Computing the current position
                    float3 current_position = input.positionOS + _StepSize * i * input.ray_directionOS;
                    
                    // If the current position left the cube boundaries, break
                    if(max(abs(current_position.x), max(abs(current_position.y), abs(current_position.z)))
                        > 0.5f + EPSILON)
                        break;

                    // Sampling the texture at the current position mapped to UV range [0,1]
                    total_density += tex3D(_DensityMap, current_position + float3(0.5f, 0.5f, 0.5f)).x
                        * _StepSize * _DensityMultiplier;
                }
 
                return float4(total_density, 1.0f);
            }
        
        ENDCG
        }

    }
}
