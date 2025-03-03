using System;
using UnityEngine;

public class RayMarchedFluid : MonoBehaviour
{
    [Header("Density Map")] 
    [SerializeField] public Vector3Int size;

    [Header("Density Slice")] 
    [SerializeField] public Material sliceMaterial;
    [SerializeField] [Range(0.0f, 1.0f)] public float sliceDepth;
    [SerializeField] public Vector2 depictedDensityRange;
    
    // Density map
    [HideInInspector] public RenderTexture densityMap;
    
    // Ray Marching Compute Shader
    private ComputeShader _rayMarchingComputeShader;
    
    // Compute Shader's variable IDs
    private readonly int _densityMapID = Shader.PropertyToID("density_map");
    
    // Slice density shader's variables IDs
    private readonly int _densityMapSliceID = Shader.PropertyToID("_DensityMap");
    private readonly int _sliceDepthID = Shader.PropertyToID("_SliceDepth");
    private readonly int _sliceMinDensityID = Shader.PropertyToID("_MinDensityValue");
    private readonly int _sliceMaxDensityID = Shader.PropertyToID("_MaxDensityValue");
    
    # region Unity Callback Functions
    
    private void Awake()
    {
        // Binding the compute shader
        _rayMarchingComputeShader = Resources.Load<ComputeShader>("RayMarching");
        
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

    private void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        Debug.Log("Test!");
    }

    private void Update()
    {
        UpdateDensitySliceVariables();
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

    #endregion
}