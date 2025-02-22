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
    private int _particlesAmount;
    
    // GPU CONSTANTS
    private const float GPUGroupSize = 256.0f;
    
    // GPU Kernels IDs
    private int _resetHistogramKernelID;
    private int _calculateHistogramKernelID;
    
    // GPU variables IDs
    private int _spatialKeysBufferID = Shader.PropertyToID("spatial_keys");
    private int _spatialIndicesBufferID = Shader.PropertyToID("spatial_indices");
    private int _spatialOffsetsBufferID = Shader.PropertyToID("spatial_offsets");

    private int _sortedKeysBufferID = Shader.PropertyToID("sorted_keys");
    private int _sortedIndicesBufferID = Shader.PropertyToID("sorted_indices");
    private int _keysHistogramBufferID = Shader.PropertyToID("keys_histogram");

    private int _targetBufferID = Shader.PropertyToID("target");
    private int _groupSumsBufferID = Shader.PropertyToID("group_sums");
    private int _targetSizeID = Shader.PropertyToID("target_size");
    
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
        
        // Setting the buffer for the resetHistogram kernel
        _spatialGridComputeShader.SetBuffer(_resetHistogramKernelID, _keysHistogramBufferID, _keysHistogramBuffer);
        _spatialGridComputeShader.SetBuffer(_resetHistogramKernelID, _spatialIndicesBufferID, SpatialIndicesBuffer);
        
        // Setting the buffer for the calculateHistogram kernel
        _spatialGridComputeShader.SetBuffer(_calculateHistogramKernelID, _spatialKeysBufferID, SpatialKeysBuffer);
        _spatialGridComputeShader.SetBuffer(_calculateHistogramKernelID, _keysHistogramBufferID, _keysHistogramBuffer);
    }
}
