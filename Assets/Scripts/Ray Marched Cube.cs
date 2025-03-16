using UnityEngine;

public class RayMarchedCube : MonoBehaviour
{
    [Header("Ray Marched Cube")] 
    public Material material;
    public float indexOfRefraction;
    
    [Header("Internal Sphere")] 
    public float sphereRadius;
    
    // Material shader's variables IDs
    private readonly int _sphereRadiusID = Shader.PropertyToID("sphere_radius");
    private readonly int _indexOfRefractionID = Shader.PropertyToID("IOR");
    
    void Update()
    {
        // Setting the variables
        material.SetFloat(_sphereRadiusID, sphereRadius);
        material.SetFloat(_indexOfRefractionID, indexOfRefraction);
    }
}
