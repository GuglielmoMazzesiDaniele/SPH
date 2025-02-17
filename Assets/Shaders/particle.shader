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
            Blend SrcAlpha OneMinusSrcAlpha
            
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma target 5.0;

            #include "UnityCG.cginc"

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };

            struct v2f
            {
                float4 vertex : SV_POSITION;
                float2 uv : TEXCOORD0;
                float4 color : COLOR;
            };

            struct particle
            {
                float2 position;
                float2 predicted_position;
                float2 velocity;
                float2 acceleration;
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
            v2f vert (appdata input, uint instanceID : SV_InstanceID)
            {
                // Initializing the output struct
                v2f output;

                // Setting the UV coordinates
                output.uv = input.uv;

                // Extracting the particle from the buffer
                particle current_particle;

                if(is_ping_active == 1)
                    current_particle = particles_buffer_pong[instanceID];
                else
                    current_particle = particles_buffer_pong[instanceID];

                // Extracting the position of the particle
                float3 particle_world_position = float3(current_particle.position, 0);
                // Translating the vertex world position by the particle position
                float3 shifted_vertex_world_position = particle_world_position + mul(unity_ObjectToWorld,
                    input.vertex * particle_radius);
                // Transforming the vertex in object space
                float3 shifted_vertex_object_position = mul(unity_WorldToObject, float4(shifted_vertex_world_position, 1));
                // Transforming the vertex in clip space, ready for rasterization
                output.vertex = UnityObjectToClipPos(shifted_vertex_object_position);

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

            // Fragment shader
            float4 frag (v2f input) : SV_Target
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
