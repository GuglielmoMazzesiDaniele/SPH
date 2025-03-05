Shader "Custom/Density Slice"
{
    Properties
    {
        // Density map input, provided by the fluid simulation
        _DensityMap ("Density Map", 3D) = "" {}
        // The slice depth to render, expressed in range [0, 1]
        _SliceDepth ("Slice Depth", Float) = 1.0
        // The minimum density value to depict
        _MinDensityValue ("Minumum Density Value", Float) = 0.0
        // The maximum density value to depict
        _MaxDensityValue ("Maximum Density Value", Float) = 5000.0
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }

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
            float _MinDensityValue;
            float _MaxDensityValue;

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

                // Initializing the color to zero
                float color = 0.0f;

                // If the density is within the accepted range, map it to range [0, 1]
                if(density >= _MinDensityValue && density <= _MaxDensityValue)
                {
                    color = (density - _MinDensityValue) / (_MaxDensityValue - _MinDensityValue);
                }

                // Convert to grayscale for visualization
                return float4(color, color, color, 1.0);
            }
            
            ENDCG
        }

    }
}
