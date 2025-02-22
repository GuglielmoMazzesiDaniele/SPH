using System.Collections.Generic;
using UnityEngine;

public class SpatialGrid
{
    // Buffers used by the fluid simulation
    public ComputeBuffer SpatialKeys;
    public ComputeBuffer SpatialIndices;
    public ComputeBuffer SpatialOffsets;
    
    // Internal auxiliary buffers used to execute exclusive prefix sum, sorting and offsets computation
    private ComputeBuffer _sortedKeys;
    private ComputeBuffer _sortedIndices;
    private ComputeBuffer _keysHistogram;
    
    // Internal pool of buffers, used to execute recursive exclusive prefix sum algorithm
    private Dictionary<int, ComputeBuffer> _buffersPool = new();
    
    // Internal reference to the spatial grid ComputeShader
    private ComputeShader _spatialGridComputeShader;
    
    // GPU variables IDs
    private int _resetHistogramAndIndicesID;

    /// <summary>
    /// Initializes the spatial grid with the given particles amount.
    /// </summary>
    /// <param name="particlesAmount">The amount of particles in the simulation.</param>
    public SpatialGrid(int particlesAmount)
    {
        // Loading the compute shader
        _spatialGridComputeShader = Resources.Load<ComputeShader>("SpatialGrid");
        
        // Initializing the spatial grid buffers
        SpatialKeys = new ComputeBuffer(particlesAmount, 3 * sizeof(float));
        SpatialIndices = new ComputeBuffer(particlesAmount, 3 * sizeof(float));
        SpatialOffsets = new ComputeBuffer(particlesAmount, 3 * sizeof(float));
    }

    /// <summary>
    /// Releases the compute buffers allocated by this class.
    /// </summary>
    public void Release()
    {
        // Releasing the public buffers
        SpatialKeys?.Release();
        SpatialIndices?.Release();
        SpatialOffsets?.Release();
        
        // Releasing the internal auxiliary buffers
        _sortedKeys?.Release();
        _sortedIndices?.Release();
        _keysHistogram?.Release();
        
        // Releasing the buffers pool
        foreach (var buffer in _buffersPool)
        {
            buffer.Value?.Release();
        }
    }
}
