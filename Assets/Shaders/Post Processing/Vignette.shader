Shader "Hidden/Vignette"
{
    Properties
    {
        // Intensity of the vignette effect
        _Intensity ("Vignette Intensity", Range(0, 1)) = 1.0
        // Smoothness of the vignette effect
        _Smoothness ("Smoothness", Range(0, 1)) = 0.0
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" "Queue"="Overlay"}

        Pass
        {
            HLSLPROGRAM
            #pragma vertex vert
            #pragma fragment frag

            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            #include "Packages/com.unity.render-pipelines.core/Runtime/Utilities/Blit.hlsl"
            
            // Variables obtained from properties
            float _Intensity;
            float _Smoothness;

            struct CustomAttributes
            {
                uint vertexID : SV_VertexID;
                UNITY_VERTEX_INPUT_INSTANCE_ID
            };

            struct CustomVaryings
            {
                float4 positionCS : SV_POSITION;
                float2 texcoord   : TEXCOORD0;
                UNITY_VERTEX_OUTPUT_STEREO
            };

            Varyings vert (CustomAttributes input)
            {
                Varyings output;
                UNITY_SETUP_INSTANCE_ID(input);
                UNITY_INITIALIZE_VERTEX_OUTPUT_STEREO(output);

                float4 pos = GetFullScreenTriangleVertexPosition(input.vertexID);
                float2 uv  = float2((input.vertexID << 1) & 2, 1.0 - (input.vertexID & 2));

                output.positionCS = pos;
                output.texcoord   = DYNAMIC_SCALING_APPLY_SCALEBIAS(uv);

                return output;
            }

            // Fragment Shader
            float4 frag(Varyings fragment) : SV_Target
            {
                // Sample the main texture (camera color)
                float3 color = SAMPLE_TEXTURE2D(_BlitTexture, sampler_LinearClamp, fragment.texcoord).rgb;
                
                // Compute vignette effect
                float2 center = float2(0.5, 0.5);
                float dist = distance(fragment.texcoord, center) * _Intensity;
                float vignette = smoothstep(1.0, _Smoothness, dist);
                
                // Apply vignette effect
                return float4(color * vignette, 1);
            }
            
            ENDHLSL
        }
    }
}
