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
            #pragma vertex Vert
            #pragma fragment frag

            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            #include "Packages/com.unity.render-pipelines.core/Runtime/Utilities/Blit.hlsl"
            
            // The resulting render of the camera up until this pass
            TEXTURE2D(_MainTex);
            float4 _MainTex_ST;
            // Sampler for the render
            SAMPLER(sampler_MainTex);

            // Variables obtained from properties
            float _Intensity;
            float _Smoothness;

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
                color *= vignette;
                
                return float4(color, 1.0);
            }


            ENDHLSL
        }
    }
}
