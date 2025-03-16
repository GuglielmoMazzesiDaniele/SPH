Shader "Custom/Floor"
{
    Properties
    {
        _FirstQuadrantColor ("First Quadrant Color", Color) = (.75, 0, 0, 1)
        _SecondQuadrantColor ("Second Quadrant Color", Color) = (0, .75, 0, 1)
        _ThirdQuadrantColor ("Third Quadrant Color", Color) = (0, 0, .75, 1)
        _FourthQuadrantColor ("Fourth Quadrant Color", Color) = (.5, .5, .5, 1)
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }

        Pass
        {
            HLSLPROGRAM
            #pragma vertex vert
            #pragma fragment frag

            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"

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
            
            // Properties
            float4 _FirstQuadrantColor;
            float4 _SecondQuadrantColor;
            float4 _ThirdQuadrantColor;
            float4 _FourthQuadrantColor;
            
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
                // Mapping the UV coordinates from [0, 1] to [-1, 1]
                float2 expanded_uv = input.uv * 2 - 1.0f;

                // Initializing the quadrant color
                float4 final_color = 0;
                
                // First quadrant 
                if(expanded_uv.x >= 0 && expanded_uv.y >= 0)
                    final_color =  _FirstQuadrantColor;
                // Second quadrant
                if(expanded_uv.x < 0 && expanded_uv.y >= 0)
                    final_color = _SecondQuadrantColor;
                // Third quadrant
                if(expanded_uv.x < 0 && expanded_uv.y < 0)
                    final_color = _ThirdQuadrantColor;
                // Fourth quadrant
                if(expanded_uv.x >= 0 && expanded_uv.y < 0)
                    final_color = _FourthQuadrantColor;

                // Chessboard pattern
                if(abs(floor(expanded_uv.x * 100) % 2) == abs(floor(expanded_uv.y * 100) % 2))
                    final_color = float4(
                        final_color.x * 0.75,
                        final_color.y * 0.75,
                        final_color.z * 0.75,
                        final_color.w);

                return final_color;
            }
            
            ENDHLSL
        }
    }
}
