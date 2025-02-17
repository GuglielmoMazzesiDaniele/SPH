using UnityEngine;

public static class Drawer
{
    
    public static void Circle(Vector3 worldPosition, float radius, Color color)
    {
        // Generating the Quad Mesh
        var quadMesh = MeshGenerator.CreateQuadMesh();

        // Creating the material of the point
        var material = new Material(Shader.Find("Custom/Particle"));

        // Set color property
        material.SetColor("_Color", color);
        
        // Create transformation matrix for the quad
        var transform = Matrix4x4.TRS(worldPosition, 
            Quaternion.identity, 
            new Vector3(radius * 2, radius * 2, 1));
        
        // Initializing the rendering parameters
        var parameters = new RenderParams(material);

        // Rendering the mesh
        Graphics.RenderMesh(parameters, quadMesh, 0, transform);
    }
}