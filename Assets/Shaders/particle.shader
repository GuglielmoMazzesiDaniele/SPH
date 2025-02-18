Shader "Custom/Particle"
{
    Properties
    {
        _Color ("Color", Color) = (1, 1, 1, 1)
        _SpeedGradientTex ("Speed Gradient Texture", 2D) = "white" {}
        _MaxVelocity ("Max Velocity", Float) = 5.0
    }
    
    SubShader
    {
        Tags { "RenderType"="Transparent" "Queue"="Overlay" }
        
        Pass
        {
            ZWrite Off
            ZTest LEqual
            Blend SrcAlpha OneMinusSrcAlpha
            
            CGPROGRAM
            #pragma vertex vert
            #pragma geometry geom
            #pragma fragment frag
            #pragma target 5.0;

            #include "UnityCG.cginc"

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };

            struct v2g
            {
                float4 clip_pos : SV_POSITION;
                float4 color : COLOR;
            };

            struct g2f
            {
                float4 clip_pos : SV_POSITION;
                float2 uv : TEXCOORD0;
                float4 color : COLOR;
            };

            struct particle
            {
                float3 position;
                float3 predicted_position;
                float3 velocity;
                float3 acceleration;
                float density;
            };
            
            // Data received from Properties
            float4 _Color;
            float _MaxVelocity;
            Texture2D _SpeedGradientTex;
            SamplerState linear_clamp_sampler;

            // CPU PROVIDED DATA
            float particle_radius;
            int is_ping_active;

            // GPU RELATED DATA
            StructuredBuffer<particle> particles_buffer_ping;
            StructuredBuffer<particle> particles_buffer_pong;
            
            // Vertex shader
            v2g vert (uint instanceID : SV_VertexID)
            {
                // Initializing the output struct
                v2g output;
                
                // Extracting the particle from the buffer
                particle current_particle = particles_buffer_ping[instanceID];

                // Extracting the particle from the active buffer
                // if(is_ping_active == 1)
                //     current_particle = particles_buffer_ping[instanceID];
                // else
                //     current_particle = particles_buffer_pong[instanceID];

                // Transforming the particle position in object space
                output.clip_pos = UnityObjectToClipPos(float4(current_particle.position, 1));
                
                // Extracting the magnitude of the particle's speed
                float speed = length(current_particle.velocity);
                // Normalizing the particle speed in range [0, 1], clamping it to 1 using saturate function
                float speed_n = saturate(speed / _MaxVelocity);
                // Computing the vertex color by sampling from the gradient texture generate by the CPU, using a linear
                // clamp sampling algorithm. (apparently Unity recognizes keywords in variables name and automatically
                // generate sampling algorithms, see https://docs.unity3d.com/Manual/SL-SamplerStates.html)
                output.color = _SpeedGradientTex.SampleLevel(linear_clamp_sampler, float2(speed_n, 0.5), 0);
                
                return output;
            }

            // Geometry shader
            [maxvertexcount(6)]
            void geom(point v2g input [1], inout TriangleStream<g2f> triangle_stream)
            {
                // Extracting the world position from the struct provided by vertex shader
                float4 clip_pos = input[0].clip_pos;

                // Initializing the offset
                float2 offset = particle_radius * float2(1, _ScreenParams.x / _ScreenParams.y);
                
                // Initializing the quad
                g2f quad[4];

                // Computing the bottom left vertex
                quad[0].uv = float2(0, 0);
                quad[0].color = input[0].color;
                quad[0].clip_pos = clip_pos + float4(-offset.x, -offset.y, 0, 0);

                // Computing the bottom right vertex
                quad[1].uv = float2(1, 0);
                quad[1].color = input[0].color;
                quad[1].clip_pos = clip_pos + float4( offset.x, -offset.y, 0, 0);

                // Computing the top left vertex
                quad[2].uv = float2(0, 1);
                quad[2].color = input[0].color;
                quad[2].clip_pos = clip_pos + float4(-offset.x,  offset.y, 0, 0);

                // Computing the top right vertex
                quad[3].uv = float2(1, 1);
                quad[3].color = input[0].color;
                quad[3].clip_pos = clip_pos + float4( offset.x,  offset.y, 0, 0);

                // Appending the bottom triangle
                triangle_stream.Append(quad[0]);
                triangle_stream.Append(quad[1]);
                triangle_stream.Append(quad[2]);

                // Appending the top triangle
                triangle_stream.Append(quad[2]);
                triangle_stream.Append(quad[3]);
                triangle_stream.Append(quad[1]);
            }
            
            // Fragment shader
            float4 frag (g2f input) : SV_Target
            {
                // Mapping UV space from [0, 1] to [-1, 1] for a smoother radius computation
                float2 centre_offset = (input.uv - 0.5) * 2;

                // Computing the squared distance of the UV coordinates
                float squared_distance = dot(centre_offset, centre_offset);

                // Computing the delta of the UV coordinates with respects to neighbours
                float delta = fwidth(squared_distance);

                // Computing the alpha value of the pixel using a smoothstep function
                float alpha = 1 - smoothstep(1 - delta, 1 + delta, squared_distance);

                // Returning the pixel value
                return float4(input.color.rgb, alpha);
            }
            ENDCG
        }
    }
}
