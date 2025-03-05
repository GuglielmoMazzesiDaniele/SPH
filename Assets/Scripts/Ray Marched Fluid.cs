using System;
using UnityEngine;
using UnityEngine.Serialization;

public class RayMarchedFluid : MonoBehaviour
{
    [Header("Density Map")] 
    [SerializeField] public Vector3Int size;

    [Header("Density Slice")] 
    [SerializeField] public Material sliceMaterial;
    [SerializeField] [Range(0.0f, 1.0f)] public float sliceDepth;
    [SerializeField] public Vector2 depictedDensityRange;

    [Header("Ray Marching Fluid")] 
    [SerializeField] public Material rayMarchingMaterial;
    [SerializeField] public int stepsAmount;
    [SerializeField] public float densityMultiplier;
    
    // Density map
    [HideInInspector] public RenderTexture densityMap;
    
    // Slice density shader's variables IDs
    private readonly int _densityMapSliceID = Shader.PropertyToID("_DensityMap");
    private readonly int _sliceDepthID = Shader.PropertyToID("_SliceDepth");
    private readonly int _sliceMinDensityID = Shader.PropertyToID("_MinDensityValue");
    private readonly int _sliceMaxDensityID = Shader.PropertyToID("_MaxDensityValue");
    
    // Ray Marching fluid shader's variables IDs
    private readonly int _stepsAmountID = Shader.PropertyToID("_StepsAmount");
    private readonly int _densityMultiplierID = Shader.PropertyToID("_DensityMultiplier");
    private readonly int _stepSizeID = Shader.PropertyToID("step_size");
    
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
        
        // Assigning the density map to the material
        sliceMaterial.SetTexture(_densityMapSliceID, densityMap);
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
        // Assigning the slice depth to the material
        sliceMaterial.SetFloat(_sliceDepthID, sliceDepth);
        
        // Setting the range of densities to depict
        sliceMaterial.SetFloat(_sliceMinDensityID, depictedDensityRange.x);
        sliceMaterial.SetFloat(_sliceMaxDensityID, depictedDensityRange.y);
    }
    
    private void UpdateRayMarchingVariables()
    {
        // Assigning the density map
        rayMarchingMaterial.SetTexture(_densityMapSliceID, densityMap);
        
        // Setting the step size
        rayMarchingMaterial.SetInt(_stepsAmountID, stepsAmount);
        rayMarchingMaterial.SetFloat(_densityMultiplierID, densityMultiplier);
        rayMarchingMaterial.SetFloat(_stepSizeID, 1.0f / stepsAmount);
        
        // Setting the range of densities to depict
        rayMarchingMaterial.SetFloat(_sliceMinDensityID, depictedDensityRange.x);
        rayMarchingMaterial.SetFloat(_sliceMaxDensityID, depictedDensityRange.y);
    }

    #endregion
}