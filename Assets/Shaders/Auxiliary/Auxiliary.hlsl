#include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"

// STRUCTS
struct Ray
{
    float3 origin;
    float3 direction;
};

// GLOBAL VARIABLES
float4x4 floor_world_to_object;



// FUNCTIONS

// Given a pair of UV coordinates, generate a random value
float rand(float2 uv)
{
    return frac(sin(dot(uv, float2(12.9898, 78.233))) * 43758.5453);
}

// Given a UV coordinate belonging to floor quad, return its corresponding color
float4 floor_uv_to_color(float2 uv)
{
    // Mapping the UV coordinates from [0, 1] to [-1, 1]
    float2 expanded_uv = uv * 2 - 1.0f;

    // Initializing the quadrant color
    float3 final_color = 0;
                
    // First quadrant 
    if(expanded_uv.x >= 0 && expanded_uv.y >= 0)
        final_color = float3(0.941, 0.627, 0.294);
    // Second quadrant
    if(expanded_uv.x < 0 && expanded_uv.y >= 0)
        final_color = float3(0.694, 0.761, 0.620);
    // Third quadrant
    if(expanded_uv.x < 0 && expanded_uv.y < 0)
        final_color = float3(0.980, 0.855, 0.478);
    // Fourth quadrant
    if(expanded_uv.x >= 0 && expanded_uv.y < 0)
        final_color =  float3(0.988, 0.906, 0.784);
    
    // Chessboard pattern
    if(abs(floor(expanded_uv.x * 100) % 2) == abs(floor(expanded_uv.y * 100) % 2))
    {
        // Generating a random darkening value between 0.7 and 1
        float random_value = rand(floor(uv * 100));

        // Mapping the random value from (0, 1] to [0.7, 1]
        random_value = 0.5 + 0.5 * random_value;
        
        // Darkening the current cell
        final_color *= random_value;
    }

    return float4(final_color, 1);
}

// Given a ray in WS and the floor inverse model matrix,
// returns the UV coordinates if its intersection with the floor quad.
// If the intersection is not available, returns (-1, -1)
float2 ray_floor_intersection(Ray ray, float4x4 floor_world_to_object)
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

// Given a position in OS within a unit cube, returns the closest face normal
float3 cube_closest_face_normal(float3 position)
{
    // Computing the absolute value of each axis to find the most extended one
    float absolute_x = abs(position.x);
    float absolute_y = abs(position.y);
    float absolute_z = abs(position.z);
                
    // Case in which X axis of the position is the largest
    if (absolute_x >= absolute_y && absolute_x >= absolute_z)
    {
        // Returning the normal of the closest face perpendicular to the X axis
        return float3(sign(position.x), 0, 0);
    }
                
    // Case in which Y axis of the position is the largest
    if (absolute_y >= absolute_x && absolute_y >= absolute_z)
    {
        // Returning the normal of the closest face perpendicular to the Y axis
        return float3(0, sign(position.y), 0);
    }

    // Returning the normal of the closest face perpendicular to the Z axis
    return float3(0, 0, sign(position.z));
}

// Given a ray in OS, sample light received from the environment
float3 sample_environment(Ray ray)
{
    // Converting the Ray in WS
    ray.origin = TransformObjectToWorld(ray.origin);
    ray.direction = normalize(TransformObjectToWorldDir(ray.direction));

    // Trying to sample the floor
    float2 floor_uv = ray_floor_intersection(ray, floor_world_to_object);

    // If the ray intersect the floor, returns it color
    if(all(floor_uv != -1))
        return floor_uv_to_color(floor_uv);
                
    // Sampling the sky
    return float3(1, 1, 1);
}