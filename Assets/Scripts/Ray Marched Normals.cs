using UnityEngine;

public class RayMarchedNormals : MonoBehaviour
{
    [Header("Ray Marched Normals")] 
    public Material material;
    public int stepsAmount;
    [Range(1, 5)] public int surfaceCollisions;
    public float refractionIndex;
    public float densityMultiplier;
    public float surfaceDensityThreshold;

    // Material shader's variables IDs
    private readonly int _stepsAmountID = Shader.PropertyToID("steps_amount");
    private readonly int _maxSurfaceCollisionsID = Shader.PropertyToID("max_surface_collisions");
    private readonly int _refractionIndexID = Shader.PropertyToID("refraction_index");
    private readonly int _densityMultiplierID = Shader.PropertyToID("density_multiplier");
    private readonly int _surfaceDensityThresholdID = Shader.PropertyToID("surface_density_threshold");
    
    // Update is called once per frame
    void Update()
    {
        // Setting the variables
        material.SetInt(_stepsAmountID, stepsAmount);
        material.SetInt(_maxSurfaceCollisionsID, surfaceCollisions);
        material.SetFloat(_densityMultiplierID, densityMultiplier);
        material.SetFloat(_refractionIndexID, refractionIndex);
        material.SetFloat(_surfaceDensityThresholdID, surfaceDensityThreshold);
    }
}
