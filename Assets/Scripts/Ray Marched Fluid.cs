using UnityEngine;

public class RayMarchedFluid : MonoBehaviour
{
    [Header("Ray Marched Fluid")] 
    public Material material;
    public int stepsAmount;
    [Range(1, 3)] public int maxSurfaceCollisions;
    public float refractionIndex;
    public float densityMultiplier;
    public Vector3 scatteringCoefficients;
    
    // Material shader's variables IDs
    private readonly int _stepsAmountID = Shader.PropertyToID("steps_amount");
    private readonly int _maxSurfaceCollisionsID = Shader.PropertyToID("max_surface_collisions");
    private readonly int _refractionIndexID = Shader.PropertyToID("refraction_index");
    private readonly int _densityMultiplierID = Shader.PropertyToID("density_multiplier");
    private readonly int _scatteringCoefficientsID = Shader.PropertyToID("scattering_coefficients");

    private void Update()
    {
        // Setting the variables
        material.SetInt(_stepsAmountID, stepsAmount);
        material.SetInt(_maxSurfaceCollisionsID, maxSurfaceCollisions);
        material.SetFloat(_densityMultiplierID, densityMultiplier);
        material.SetFloat(_refractionIndexID, refractionIndex);
        material.SetVector(_scatteringCoefficientsID, scatteringCoefficients);
    }
}
