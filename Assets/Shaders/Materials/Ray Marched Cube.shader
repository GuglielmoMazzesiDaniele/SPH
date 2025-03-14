Shader "Custom/Ray Marched Cube"
{
    Properties
    {
        _Radius ("Radius", Float) = 0.5
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        Pass
        {
            HLSLPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma target 5.0

            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            
            // Allowed floating point inaccuracy
            #define EPSILON 1e-5

            struct Attributes
            {
                float4 positionOS : POSITION;
            };
            
            struct Varyings
            {
                float4 positionCS : SV_POSITION;
                float3 positionOS: TEXCOORD1;
                float3 view_directionWS: TEXCOORD2;
            };

            Varyings vert (Attributes input)
            {
                Varyings output;

                // Computing the output
                output.positionCS = TransformObjectToHClip(input.positionOS);
                output.positionOS = input.positionOS;

                // Computing the view direction in world space
                output.view_directionWS = - GetWorldSpaceViewDir(TransformObjectToWorld(input.positionOS));
                
                return output;
            }

            float _Radius;
            
            float4 frag (Varyings input) : SV_Target
            {
                // Computing the view direction
                float3 view_directionOS = TransformWorldToObjectDir(input.view_directionWS);

                // Initializing the color
                float4 color = 0;

                // Initializing the position
                float3 current_position = input.positionOS;

                // Centre of the sphere
                float3 centre = 0;

                for(int i = 0; i < 250; i++)
                {
                    if(distance(current_position, centre) <= _Radius)
                    {
                        color = float4(1, 0, 0, 1);
                        break;
                    }

                    current_position += view_directionOS * (1.0 / 100);
                }
                
                return  color;
            }
            ENDHLSL
        }
    }
}
