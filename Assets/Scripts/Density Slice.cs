using UnityEngine;
using UnityEngine.Serialization;

public class DensitySlice : MonoBehaviour
{
    [Header("Density Slice")] 
    public Material material;
    [Range(0.0f, 1.0f)] public float sliceDepth;
    public Vector2 depictedDensityRange;
    
    // Material shader's variables IDs
    private readonly int _sliceDepthID = Shader.PropertyToID("_SliceDepth");
    private readonly int _sliceMinDensityID = Shader.PropertyToID("_MinDensityValue");
    private readonly int _sliceMaxDensityID = Shader.PropertyToID("_MaxDensityValue");
    
    void Update()
    {
        // Setting the variables
        material.SetFloat(_sliceDepthID, sliceDepth);
        material.SetFloat(_sliceMinDensityID, depictedDensityRange.x);
        material.SetFloat(_sliceMaxDensityID, depictedDensityRange.y);
    }
}
