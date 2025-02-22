using System.Collections.Generic;
using UnityEngine;

public class SpatialGrid
{
    // Buffers used by the fluid simulation
    public ComputeBuffer SpatialKeysBuffer;
    public ComputeBuffer SpatialIndicesBuffer;
    public ComputeBuffer SpatialOffsetsBuffer;
    
    // Internal auxiliary buffers used to execute exclusive prefix sum, sorting and offsets computation
    private ComputeBuffer _sortedKeysBuffer;
    private ComputeBuffer _sortedIndicesBuffer;
    public ComputeBuffer _keysHistogramBuffer;
    
    // Internal pool of buffers, used to execute recursive exclusive prefix sum algorithm
    private Dictionary<int, ComputeBuffer> _buffersPool = new();
    
    // Internal reference to the spatial grid ComputeShader
    private ComputeShader _spatialGridComputeShader;
    
    // Amount of particles currently active in the simulation
    private readonly int _particlesAmount;
    
    // GPU CONSTANTS
    private const float GPUGroupSize = 256.0f;
    
    // GPU Kernels IDs
    private int _resetHistogramKernelID;
    private int _calculateHistogramKernelID;
    private int _espScanKernelID;
    private int _espCombineKernelID;
    
    // GPU variables IDs
    private readonly int _spatialKeysBufferID = Shader.PropertyToID("spatial_keys");
    private readonly int _spatialIndicesBufferID = Shader.PropertyToID("spatial_indices");
    private int _spatialOffsetsBufferID = Shader.PropertyToID("spatial_offsets");

    private int _sortedKeysBufferID = Shader.PropertyToID("sorted_keys");
    private int _sortedIndicesBufferID = Shader.PropertyToID("sorted_indices");
    private readonly int _keysHistogramBufferID = Shader.PropertyToID("keys_histogram");

    private readonly int _espTargetBufferID = Shader.PropertyToID("esp_target");
    private readonly int _espGroupSumsBufferID = Shader.PropertyToID("esp_groups_sums");
    private readonly int _espTargetSizeID = Shader.PropertyToID("esp_target_size");
    
    /// <summary>
    /// Initializes the spatial grid with the given particles amount.
    /// </summary>
    /// <param name="particlesAmount">The amount of particles in the simulation.</param>
    public SpatialGrid(int particlesAmount)
    {
        // Loading the compute shader
        _spatialGridComputeShader = Resources.Load<ComputeShader>("SpatialGrid");
        
        // Storing the particles amount
        _particlesAmount = particlesAmount;
        
        // Initializing the buffers 
        InitializeBuffers();
        
        // Initializing the compute shaders variables
        InitializeComputeShader();
    }
    
    /// <summary>
    /// Updates the spatial lookup buffers by using the data contained within the SpatialKeys buffer.
    /// </summary>
    public void UpdateSpatialLookup()
    {
        // Resetting the histogram and indices buffers
        _spatialGridComputeShader.Dispatch(_resetHistogramKernelID, 
            Mathf.CeilToInt(_particlesAmount / GPUGroupSize), 1, 1);
        
        // Computing the histogram of the keys value
        _spatialGridComputeShader.Dispatch(_calculateHistogramKernelID,
            Mathf.CeilToInt(_particlesAmount / GPUGroupSize), 1, 1);
        
        // Applying ESP to keys histogram buffer
        ExclusivePrefixSum(_keysHistogramBuffer);
    }
    
    /// <summary>
    /// Releases the compute buffers allocated by this class.
    /// </summary>
    public void Release()
    {
        // Releasing the public buffers
        SpatialKeysBuffer?.Release();
        SpatialIndicesBuffer?.Release();
        SpatialOffsetsBuffer?.Release();
        
        // Releasing the internal auxiliary buffers
        _sortedKeysBuffer?.Release();
        _sortedIndicesBuffer?.Release();
        _keysHistogramBuffer?.Release();
        
        // Releasing the buffers pool
        foreach (var buffer in _buffersPool)
        {
            buffer.Value?.Release();
        }
    }

    /// <summary>
    /// Given a target buffer, executes the exclusive prefix sum on it.
    /// CAREFUL: this function recursively computes the effects of GPU groups to sequentials groups.
    /// </summary>
    /// <param name="target"></param>
    private void ExclusivePrefixSum(ComputeBuffer target)
    {
        // Computing the number of groups needed to handle the target buffer
        var requiredGroupsAmount = Mathf.CeilToInt(target.count / 2.0f / 256.0f);
        
        // Trying to obtain a buffer from the pool having the required size
        if (!_buffersPool.TryGetValue(requiredGroupsAmount, out var espGroupsSums))
        {
            // Initializing the esp groups sums buffer
            espGroupsSums = new ComputeBuffer(requiredGroupsAmount, sizeof(uint));
            // Adding the buffer to the pool
            _buffersPool.Add(requiredGroupsAmount, espGroupsSums);
        }
        
        // Setting the buffers to the esp scan kernel
        _spatialGridComputeShader.SetBuffer(_espScanKernelID, _espTargetBufferID, target);
        _spatialGridComputeShader.SetBuffer(_espScanKernelID, _espGroupSumsBufferID, espGroupsSums);
        
        // Setting the current esp target size
        _spatialGridComputeShader.SetInt(_espTargetSizeID, target.count);
        
        // Dispatching the esp scan kernel
        _spatialGridComputeShader.Dispatch(_espScanKernelID, requiredGroupsAmount, 1, 1);
        
        // If more than one group was required to execute the esp algorithm, the total sum of each group must be
        // increased by the total sum of all previous groups. To do so, esp is applied to the group sum buffer as well.
        if (requiredGroupsAmount > 1)
        {
            // Recursively calculate esp also on the buffer containing the total sum of each group
            ExclusivePrefixSum(espGroupsSums);
            
            // Setting the buffers to the esp combine kernel
            _spatialGridComputeShader.SetBuffer(_espCombineKernelID, _espTargetBufferID, target);
            _spatialGridComputeShader.SetBuffer(_espCombineKernelID, _espGroupSumsBufferID, espGroupsSums);
            
            // Setting the current esp target size
            _spatialGridComputeShader.SetInt(_espTargetSizeID, target.count);
            
            // Dispatching the esp combine kernel
            _spatialGridComputeShader.Dispatch(_espCombineKernelID, requiredGroupsAmount, 1, 1);
        }
    }

    /// <summary>
    /// Initializes the compute buffers
    /// </summary>
    private void InitializeBuffers()
    {
        // Initializing the spatial grid buffers
        SpatialKeysBuffer = new ComputeBuffer(_particlesAmount, sizeof(uint));
        SpatialIndicesBuffer = new ComputeBuffer(_particlesAmount, sizeof(uint));
        SpatialOffsetsBuffer = new ComputeBuffer(_particlesAmount, sizeof(uint));
        
        // Initializing the auxiliary buffers
        _sortedKeysBuffer = new ComputeBuffer(_particlesAmount, sizeof(uint));
        _sortedIndicesBuffer = new ComputeBuffer(_particlesAmount, sizeof(uint));
        _keysHistogramBuffer = new ComputeBuffer(_particlesAmount, sizeof(uint));
    }
    
    /// <summary>
    /// Initializes the compute shader variables.
    /// </summary>
    private void InitializeComputeShader()
    {
        // Setting the compute shader ints
        _spatialGridComputeShader.SetInt("particles_amount", _particlesAmount);
        
        // Setting the compute shaders kernel IDs
        _resetHistogramKernelID = _spatialGridComputeShader.FindKernel("reset_histogram");
        _calculateHistogramKernelID = _spatialGridComputeShader.FindKernel("calculate_histogram");
        _espScanKernelID = _spatialGridComputeShader.FindKernel("esp_scan");
        _espCombineKernelID = _spatialGridComputeShader.FindKernel("esp_combine");
        
        // Setting the buffers for the resetHistogram kernel
        _spatialGridComputeShader.SetBuffer(_resetHistogramKernelID, _keysHistogramBufferID, _keysHistogramBuffer);
        _spatialGridComputeShader.SetBuffer(_resetHistogramKernelID, _spatialIndicesBufferID, SpatialIndicesBuffer);
        
        // Setting the buffers for the calculateHistogram kernel
        _spatialGridComputeShader.SetBuffer(_calculateHistogramKernelID, _spatialKeysBufferID, SpatialKeysBuffer);
        _spatialGridComputeShader.SetBuffer(_calculateHistogramKernelID, _keysHistogramBufferID, _keysHistogramBuffer);
    }
}
