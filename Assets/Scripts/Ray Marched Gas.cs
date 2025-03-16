using UnityEngine;

public class RayMarchedGas : MonoBehaviour
{
    [Header("Ray Marching Gas")] 
    public Material material;
    public int stepsAmount;
    public int internalStepsAmount;
    public Vector3 sunDirection;
    public float densityMultiplier;
    public Vector3 scatteringCoefficients;
    
    // Material shader's variables IDs
    private readonly int _stepsAmountID = Shader.PropertyToID("steps_amount");
    private readonly int _internalStepsAmountID = Shader.PropertyToID("internal_steps_amount");
    private readonly int _sunDirectionID = Shader.PropertyToID("sun_direction");
    private readonly int _densityMultiplierID = Shader.PropertyToID("density_multiplier");
    private readonly int _scatteringCoefficientsID = Shader.PropertyToID("scattering_coefficients");
    
    private void Update()
    {
        // Setting the variables
        material.SetInt(_stepsAmountID, stepsAmount);
        material.SetInt(_internalStepsAmountID, internalStepsAmount);
        material.SetVector(_sunDirectionID, sunDirection.normalized);
        material.SetFloat(_densityMultiplierID, densityMultiplier);
        material.SetVector(_scatteringCoefficientsID, scatteringCoefficients);
    }
}