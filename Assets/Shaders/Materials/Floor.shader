Shader "Custom/Floor"
{
    SubShader
    {
        Tags { "RenderType"="Opaque" }

        Pass
        {
            HLSLPROGRAM
            #pragma vertex vert
            #pragma fragment frag

            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            #include "../Auxiliary/Auxiliary.hlsl"

            struct Attributes
            {
                float4 positionOS : POSITION;
                float2 uv : TEXCOORD0;
            };
            
            struct Varyings
            {
                float4 positionCS : SV_POSITION;
                float2 uv : TEXCOORD0;
                float3 positionWS: TEXCOORD1;
            };
            
            // Vertex shader
            Varyings vert (Attributes input)
            {
                // Initializing the output
                Varyings output;

                // Computing the vertex position in clip space
                output.positionCS = TransformObjectToHClip(input.positionOS.xyz);

                // Storing the position in object space
                output.positionWS = TransformObjectToWorld(input.positionOS);

                // Storing the UV coordinates
                output.uv = input.uv;
                
                return output;
            }

            // Fragment shader
            float4 frag (Varyings input) : SV_Target
            {
                return floor_uv_to_color(input.uv);
            }
            
            ENDHLSL
        }
    }
}
