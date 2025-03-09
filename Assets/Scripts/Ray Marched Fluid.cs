using System;
using UnityEngine;
using UnityEngine.Serialization;

public class RayMarchedFluid : MonoBehaviour
{
    [Header("Density Map")] 
    [SerializeField] public Vector3Int size;

    [Header("Density Slice")] 
    public Material sliceMaterial;
    [Range(0.0f, 1.0f)] public float sliceDepth;
    public Vector2 depictedDensityRange;

    [Header("Ray Marching Fluid")] 
    public Material rayMarchingMaterial;
    public int stepsAmount;
    public int internalStepsAmount;
    public Vector3 sunDirection;
    public float densityMultiplier;
    public Vector3 scatteringCoefficients;
    
    // Density map
    [HideInInspector] public RenderTexture densityMap;
    
    // Slice density shader's variables IDs
    private readonly int _densityMapSliceID = Shader.PropertyToID("_DensityMap");
    private readonly int _sliceDepthID = Shader.PropertyToID("_SliceDepth");
    private readonly int _sliceMinDensityID = Shader.PropertyToID("_MinDensityValue");
    private readonly int _sliceMaxDensityID = Shader.PropertyToID("_MaxDensityValue");
    
    // Ray Marching fluid shader's variables IDs
    private readonly int _stepsAmountID = Shader.PropertyToID("steps_amount");
    private readonly int _internalStepsAmountID = Shader.PropertyToID("internal_steps_amount");
    private readonly int _sunDirectionID = Shader.PropertyToID("sun_direction");
    private readonly int _densityMultiplierID = Shader.PropertyToID("density_multiplier");
    private readonly int _scatteringCoefficientsID = Shader.PropertyToID("scattering_coefficients");
    
    # region Unity Callback Functions
    
    private void Awake()
    {
        // Initializing the density map
        densityMap = new RenderTexture(size.x, size.y, 0, RenderTextureFormat.RFloat)
        {
            // Signaling to Unity that this is a 3D texture
            dimension = UnityEngine.Rendering.TextureDimension.Tex3D,
            // Setting the texture depth 
            volumeDepth = size.z,
            // Allowing random write on the texture (done by the Compute Shader)
            enableRandomWrite = true,
            // Setting the texture filter
            filterMode = FilterMode.Bilinear,
            // Setting the wrap mode
            wrapMode = TextureWrapMode.Clamp
        };
        
        // Creating the density map
        densityMap.Create();
        
        // Assigning the density map to the materials
        sliceMaterial.SetTexture(_densityMapSliceID, densityMap);
        rayMarchingMaterial.SetTexture(_densityMapSliceID, densityMap);
    }

    private void Update()
    {
        // Updating the density slice variables
        UpdateDensitySliceVariables();
        
        // Updating the ray marching fluid variable
        UpdateRayMarchingVariables();
    }

    /// <summary>
    /// Updates the variable used by the density slice material
    /// </summary>
    private void UpdateDensitySliceVariables()
    {
        // Setting the variables
        sliceMaterial.SetFloat(_sliceDepthID, sliceDepth);
        sliceMaterial.SetFloat(_sliceMinDensityID, depictedDensityRange.x);
        sliceMaterial.SetFloat(_sliceMaxDensityID, depictedDensityRange.y);
    }
    
    /// <summary>
    /// Updates the variable used by the ray marching material
    /// </summary>
    private void UpdateRayMarchingVariables()
    {
        // Setting the variables
        rayMarchingMaterial.SetInt(_stepsAmountID, stepsAmount);
        rayMarchingMaterial.SetInt(_internalStepsAmountID, internalStepsAmount);
        rayMarchingMaterial.SetVector(_sunDirectionID, sunDirection.normalized);
        rayMarchingMaterial.SetFloat(_densityMultiplierID, densityMultiplier);
        rayMarchingMaterial.SetVector(_scatteringCoefficientsID, scatteringCoefficients);
    }

    #endregion
}