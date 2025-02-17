using UnityEngine;

public static class MeshGenerator
{
    /// <summary>
    /// Generates a simple 2D quad mesh
    /// </summary>
    /// <returns>The quad mesh</returns>
    public static Mesh CreateQuadMesh()
    {
        // Initializing the mesh
        var quadMesh = new Mesh();

        // Initializing the coordinates of the vertices in object space
        Vector3 [] vertices = {
            new (-0.5f, -0.5f, 0),
            new ( 0.5f, -0.5f, 0),
            new (-0.5f,  0.5f, 0),
            new ( 0.5f,  0.5f, 0)
        };

        // Array of triangles, expressed as clockwise indexes of vertices
        int [] triangles = { 0, 2, 1, 1, 2, 3 };
        // Array of UV coordinates
        Vector2 [] uvs =
        {
            new (0, 0), 
            new (1, 0), 
            new (0, 1), 
            new (1, 1)
        };

        // Initializing the mesh fields
        quadMesh.vertices = vertices;
        quadMesh.triangles = triangles;
        quadMesh.uv = uvs;

        // Returning the mesh
        return quadMesh;
    }
}