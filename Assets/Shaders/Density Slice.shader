Shader "Custom/Density Slice"
{
    Properties
    {
        // Density map input, provided by the fluid simulation
        _DensityMap ("Density Map", 3D) = "" {}
        // The slice depth to render, expressed in range [0, 1]
        _SliceDepth ("Slice Depth", Float) = 1.0
    }
    SubShader
    {
        Tags { "Queue" = "Overlay" "RenderType"="Opaque" }

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma target 5.0

            #include "UnityCG.cginc"

            struct appdata_t
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };

            struct v2f
            {
                float4 position : SV_POSITION;
                float2 uv : TEXCOORD0;
            };

            // Texture and Uniforms
            sampler3D _DensityMap;
            float _SliceDepth;

            // Vertex Shader
            v2f vert (appdata_t v)
            {
                v2f output;
                output.position = UnityObjectToClipPos(v.vertex);
                output.uv = v.uv;
                return output;
            }

            // Fragment Shader (Extracting the Slice)
            fixed4 frag (v2f i) : SV_Target
            {
                // Converting UV coordinates to UVW using the provided depth
                float3 uvw = float3(i.uv, _SliceDepth);
                
                // Sampling the 3D texture at the computed UVW coordinate
                float density = tex3D(_DensityMap, uvw).r;

                density /= 650.0f;

                // Convert to grayscale for visualization
                return float4(density, density, density, 1.0);
            }
            
            ENDCG
        }

    }
}
