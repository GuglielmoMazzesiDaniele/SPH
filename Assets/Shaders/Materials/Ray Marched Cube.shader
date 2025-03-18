Shader "Custom/Ray Marched Cube"
{
    SubShader
    {
        Tags {"RenderType"="Opaque" }
        Pass
        {
            HLSLPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma target 5.0

            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            #include "../Auxiliary/Auxiliary.hlsl"
            
            // Allowed floating point inaccuracy
            #define EPSILON 1e-5
            #define STEPS_AMOUNT 1000
            #define AIR_IOR 1.0003

            // Properties
            float sphere_radius;
            float IOR;

            // Variables
            // float4x4 floor_world_to_object;

            // Structs
            
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

            float4 frag (Varyings input) : SV_Target
            {
                // Computing the view direction
                float3 ray_direction = TransformWorldToObjectDir(input.view_directionWS);

                // Initializing the color
                float3 final_color = 0;

                // Initializing the position
                float3 current_position = input.positionOS;

                // Computing the closest face normal
                float3 normal = cube_closest_face_normal(current_position);

                // Verifying if the direction and normal are correctly aligned
                float3 correct_normal = dot(ray_direction, normal) < 0 ?
                    normal : -normal;

                // Verifying the index of refraction between the mediums
                float correct_IOR = dot(ray_direction, normal) < 0 ?
                    AIR_IOR / IOR : IOR / AIR_IOR;
                
                // Computing the refraction direction
                float3 refracted_direction = refract(ray_direction, correct_normal, correct_IOR);

                // If refraction is not possible return white
                if(all(refracted_direction == 0))
                    return float4(1,1,1,1);

                // Assigning the new directon
                ray_direction = refracted_direction;

                [loop]
                for(int i = 0; i < STEPS_AMOUNT; i++)
                {
                    // Verifying intersection with the internal sphere
                    if(distance(current_position, float3(0, 0, 0)) <= sphere_radius)
                    {
                        // Computing the normal at the current sphere surface point
                        normal = normalize(current_position);
                    
                        // Applying a basic lambertian shading to the sphere
                        final_color = float3(1, 0, 0) * smoothstep(0, 0.75, dot(ray_direction, -normal));
                        break;
                    }

                    // If I reached the limit of the cube, apply refraction leaving the cube and sample the sky
                    if(any(abs(current_position) >= 0.5f + EPSILON))
                    {
                        // Computing the closest face normal
                        normal = cube_closest_face_normal(current_position);

                        // Verifying if the direction and normal are correctly aligned
                        correct_normal = dot(ray_direction, normal) < 0 ?
                            normal : -normal;

                        // Verifying the index of refraction between the mediums
                        correct_IOR = dot(ray_direction, normal) < 0 ?
                            AIR_IOR / IOR : IOR / AIR_IOR;
                        
                        // Computing the refraction direction
                         refracted_direction = refract(ray_direction, correct_normal, correct_IOR);

                        // If refraction is not possible, break
                        if(all(refracted_direction == 0))
                            break;

                        // Initializing the ray
                        Ray exiting_ray;
                        exiting_ray.origin = current_position;
                        exiting_ray.direction = refracted_direction; 

                        // Printing the UV coordinates
                        final_color = sample_environment(exiting_ray); 

                        // Exiting the raymarching loop
                        break;
                    }

                    current_position += ray_direction * (sqrt(3) / STEPS_AMOUNT);
                }
                
                return float4(final_color, 1);
            }
            ENDHLSL
        }
    }
}
