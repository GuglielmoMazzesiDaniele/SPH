using UnityEngine;
using UnityEngine.Serialization;

public class RayMarchingFluid : MonoBehaviour
{
    [Header("Density Map")] 
    [SerializeField] public Vector3Int size;

    [Header("Density Slice")] 
    [SerializeField] public Material sliceMaterial;
    [SerializeField] [Range(0.0f, 1.0f)] public float sliceDepth;
    
    // Density map
    [HideInInspector] public RenderTexture densityMap;
    
    // Ray Marching Compute Shader
    private ComputeShader _rayMarchingComputeShader;
    
    // Compute Shader's variable IDs
    private readonly int _densityMapID = Shader.PropertyToID("density_map");
    
    // Slice density shader's variables IDs
    private readonly int _densityMapSliceID = Shader.PropertyToID("_DensityMap");
    private readonly int _sliceDepthID = Shader.PropertyToID("_SliceDepth");
    
    # region Unity Callback Functions

    private void Start()
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
    }

    private void Update()
    {
        // // Create a readable CPU copy of the 3D texture
        // var debugTexture = new Texture3D(densityMap.width, densityMap.height, densityMap.volumeDepth, TextureFormat.RFloat, false);
        //
        // Graphics.CopyTexture(densityMap, debugTexture);
        //
        // // Read and print some values from the CPU copy
        // var value = debugTexture.GetPixel(0, 0, 0);
        // Debug.Log($"Density value at (0,0,0): {value.r}");
        
        // Assigning the density map to the material
        sliceMaterial.SetTexture(_densityMapSliceID, densityMap);
        
        // Assigning the slice depth to the material
        sliceMaterial.SetFloat(_sliceDepthID, sliceDepth);
    }

    #endregion
}