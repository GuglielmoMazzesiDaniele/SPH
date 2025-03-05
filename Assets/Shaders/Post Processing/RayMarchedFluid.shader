Shader "Custom/Ray Marching Fluid"
{
    Properties
    {
        // Density map input, provided by the fluid simulation
        _DensityMap ("Density Map", 3D) = "" {}
        // The ray marching step size
        _StepSize ("Step Size", Float) = 0.01
        // Simulation boundaries
        _SimulationBoundaries ("Simulation Boundaries", Vector) = (0, 0, 0, 0)
        // The minimum density value to depict
        _MinDensityValue ("Minumum Density Value", Float) = 150.0
        // The maximum density value to depict
        _MaxDensityValue ("Maximum Density Value", Float) = 1000.0
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" "Queue"="Overlay"}
        
        Pass
        {
            HLSLPROGRAM
            #pragma vertex Vert
            #pragma fragment frag
            #pragma target 5.0

            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            #include "Packages/com.unity.render-pipelines.core/Runtime/Utilities/Blit.hlsl"
            
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

            // Declaring the density map
            TEXTURE3D(_DensityTex);
            SAMPLER(sampler_DensityTex);

            // Vertex shader
            Varyings vert (CustomAttributes input)
            {
                Varyings output;
                UNITY_SETUP_INSTANCE_ID(input);
                UNITY_INITIALIZE_VERTEX_OUTPUT_STEREO(output);

                float4 pos = GetFullScreenTriangleVertexPosition(input.vertexID);
                float2 uv  = GetFullScreenTriangleTexCoord(input.vertexID);

                output.positionCS = pos;
                output.texcoord   = DYNAMIC_SCALING_APPLY_SCALEBIAS(uv);

                return output;
            }

            // Fragment Shader
            float4 frag(Varyings fragment) : SV_Target
            {
                float4 density = SAMPLE_TEXTURE3D(_DensityTex, sampler_DensityTex, float3(fragment.texcoord.xy, 0.5));
                
                return density;
            }
            
            ENDHLSL
        }
    }
}
