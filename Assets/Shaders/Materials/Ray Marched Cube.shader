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
            
            // Allowed floating point inaccuracy
            #define EPSILON 1e-5
            #define STEPS_AMOUNT 1000
            #define AIR_IOR 1.0003

            // Properties
            float sphere_radius;
            float IOR;

            // Variables
            float4x4 floor_world_to_object;
            
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

            struct Ray
            {
                float3 origin;
                float3 direction;
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

            // Given a position in OS, returns the closest face normal
            float3 closest_face_normal(float3 positionOS)
            {
                // Computing the absolute value of each axis to find the most extended one
                float absolute_x = abs(positionOS.x);
                float absolute_y = abs(positionOS.y);
                float absolute_z = abs(positionOS.z);
                
                // Case in which X axis of the position is the largest
                if (absolute_x >= absolute_y && absolute_x >= absolute_z)
                {
                    // Returning the normal of the closest face perpendicular to the X axis
                    return float3(sign(positionOS.x), 0, 0);
                }
                
                // Case in which Y axis of the position is the largest
                if (absolute_y >= absolute_x && absolute_y >= absolute_z)
                {
                    // Returning the normal of the closest face perpendicular to the Y axis
                    return float3(0, sign(positionOS.y), 0);
                }

                // Returning the normal of the closest face perpendicular to the Z axis
                return float3(0, 0, sign(positionOS.z));
            }

            // Given a ray in WS, returns the UV coordinates if its intersection with the floor quad.
            // If the intersection is not available, returns (-1, -1)
            float2 ray_floor_intersection(Ray ray)
            {
                // Constant
                const float3 quad_normal = float3(0, 0, -1);
                
                // Transforming the ray in quad OS
                ray.origin = mul(floor_world_to_object, float4(ray.origin, 1.0)).xyz;
                ray.direction = normalize(mul(floor_world_to_object, float4(ray.direction, 0.0)).xyz);
                
                // Computing the ray magnitude to intersect the quad plane
                float magnitude = - dot(ray.origin, quad_normal) / dot(ray.direction, quad_normal);

                // Verifying that the magnitude is valid
                if(magnitude <= 0)
                    return float2(-1, -1);

                // Computing the intersection point in quad OS
                float3 intersection = ray.origin + magnitude * ray.direction;
                
                // Verifying that the intersection point lies on the quad
                if(any(abs(intersection) > 0.5))
                    return float2(-1, -1);

                // Returning the UV coordinate corresponding intersection point
                return float2(intersection.x, intersection.y) + 0.5;
            }

            // Given a UV coordinate belonging to floor quad, return its corresponding color
            float4 floor_uv_to_color(float2 uv)
            {
                // Mapping the UV coordinates from [0, 1] to [-1, 1]
                float2 expanded_uv = uv * 2 - 1.0f;

                // Initializing the quadrant color
                float4 final_color = 0;
                
                // First quadrant 
                if(expanded_uv.x >= 0 && expanded_uv.y >= 0)
                    final_color =  float4(.75, 0, 0, 1);
                // Second quadrant
                if(expanded_uv.x < 0 && expanded_uv.y >= 0)
                    final_color = float4(0, .75, 0, 1);
                // Third quadrant
                if(expanded_uv.x < 0 && expanded_uv.y < 0)
                    final_color = float4(0, 0, .75, 1);
                // Fourth quadrant
                if(expanded_uv.x >= 0 && expanded_uv.y < 0)
                    final_color = float4(.5, .5, .5, 1);

                // Chessboard pattern
                if(abs(floor(expanded_uv.x * 100) % 2) == abs(floor(expanded_uv.y * 100) % 2))
                    final_color = float4(
                        final_color.x * 0.75,
                        final_color.y * 0.75,
                        final_color.z * 0.75,
                        final_color.w);

                return final_color;
            }

            float4 frag (Varyings input) : SV_Target
            {
                // Computing the view direction
                float3 ray_direction = TransformWorldToObjectDir(input.view_directionWS);

                // Initializing the color
                float4 final_color = float4(0, 0, 0, 0);

                // Initializing the position
                float3 current_position = input.positionOS;

                // Computing the closest face normal
                float3 normal = closest_face_normal(current_position);

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
                        final_color = float4(1, 0, 0, 1) * smoothstep(0, 0.75, dot(ray_direction, -normal));
                        break;
                    }

                    // If I reached the limit of the cube, apply refraction leaving the cube and sample the sky
                    if(any(abs(current_position) >= 0.5f + EPSILON))
                    {
                        // Computing the closest face normal
                        normal = closest_face_normal(current_position);

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
                        exiting_ray.origin = TransformObjectToWorld(current_position - 1e-2 * ray_direction);
                        exiting_ray.direction = TransformObjectToWorldDir(refracted_direction);

                        // Computing the UV coordinates of intersection with the floor quad
                        float2 floor_uv = ray_floor_intersection(exiting_ray);

                        // Verifying that the UVs are valid
                        if(all(floor_uv == -1))
                            break;
                        
                        // Printing the UV coordinates
                        final_color = floor_uv_to_color(floor_uv);

                        // Exiting the raymarching loop
                        break;
                    }

                    current_position += ray_direction * (sqrt(3) / STEPS_AMOUNT);
                }
                
                return final_color;
            }
            ENDHLSL
        }
    }
}
