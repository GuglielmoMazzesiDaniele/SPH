using System.Collections.Generic;
using System.Threading.Tasks;
using UnityEngine;
using UnityEngine.Serialization;
using Random = UnityEngine.Random;
using Vector2 = UnityEngine.Vector2;
using Vector3 = UnityEngine.Vector3;

public class FluidSimulation : MonoBehaviour
{
    [Header("Particle")]
    [SerializeField] public float particleRadius;
    [SerializeField] public float particleMass;
    [SerializeField] public float collisionDamping;
    
    [Header("Kernel")]
    [SerializeField] public float kernelRadius;
    
    [Header("Bounds")]
    [SerializeField] public Vector3 bounds;
    
    [Header("Simulation")]
    [SerializeField] public float gravity;
    [SerializeField] public float restDensity;
    [SerializeField] public float stiffness;
    [SerializeField] public float viscosity;
    [SerializeField] public uint simulationSubSteps;

    [Header("Optimization")] 
    [SerializeField] public float spatialGridCellSize;
    
    [Header("Rendering")]
    [SerializeField] public Material particleMaterial;
    [SerializeField] public Gradient particleSpeedGradient;
    [SerializeField] public int speedGradientTextureWidth;
    [SerializeField] public bool renderSimulation;
    
    #region CPU Related Fields

    public enum ParticlesInitialization
    {
        Grid,
        Random
    }
    
    private struct SpatialEntry
    {
        public int particleIndex;
        public int cellKey;
    }
    
    // Reference to the simulation compute shader
    private ComputeShader _simulationComputeShader;
    
    // Reference to the particles spawner
    private FluidSpawner _fluidSpawner;
    
    // Particles related fields
    private int _particlesAmount;
    private Vector3[] _positions;
    private Vector3[] _predictedPositions;
    private Vector3[] _velocities;
    private float[] _densities;
    
    // Bounds related fields
    private Vector3 _halfBounds;
    
    // Simulation related fields
    private Vector3 _gravityAcceleration;
    
    // Kernel related fields
    private float _sqrKernelRadius;
    private float _poly6Normalization;
    private float _spikyGradientNormalization;
    private float _viscosityLaplacianNormalization;
    
    // Optimization related fields
    private SpatialEntry[] _spatialLookup;
    private int[] _lookupStartIndex;
    private List<SpatialEntry>[] _buckets;
    
    // Array of neighbours particles
    private static readonly Vector3Int[] CellOffsets =
    {
        new (-1,  1,  1),
        new ( 0,  1,  1),
        new ( 1,  1,  1),
        new (-1,  0,  1),
        new ( 0,  0,  1),
        new ( 1,  0,  1),
        new (-1, -1,  1),
        new ( 0, -1,  1),
        new ( 1, -1,  1),
        new (-1,  1,  0),
        new ( 0,  1,  0),
        new ( 1,  1,  0),
        new (-1,  0,  0),
        new ( 0,  0,  0),
        new ( 1,  0,  0),
        new (-1, -1,  0),
        new ( 0, -1,  0),
        new ( 1, -1,  0),
        new (-1,  1, -1),
        new ( 0,  1, -1),
        new ( 1,  1, -1),
        new (-1,  0, -1),
        new ( 0,  0, -1),
        new ( 1,  0, -1),
        new (-1, -1, -1),
        new ( 0, -1, -1),
        new ( 1, -1, -1),
    };
    
    // Shaders related fields
    private readonly int _halfBoundsID = Shader.PropertyToID("half_bounds");
    private readonly int _collisionDampingID = Shader.PropertyToID("collision_damping");

    private readonly int _worldToLocalMatrixID = Shader.PropertyToID("world_to_local");
    private readonly int _localToWorldMatrixID = Shader.PropertyToID("local_to_world");

    private readonly int _deltaTimeID = Shader.PropertyToID("delta_time");
    
    private readonly int _gravityID = Shader.PropertyToID("gravity");
    
    private readonly int _particleRadiusID = Shader.PropertyToID("particle_radius");
    private readonly int _kernelRadiusID = Shader.PropertyToID("kernel_radius");
    private readonly int _sqrKernelRadiusID = Shader.PropertyToID("sqr_kernel_radius");
    
    private readonly int _poly6NormalizationID = Shader.PropertyToID("poly6_normalization");
    private readonly int _spikyGradientFirstTermID = Shader.PropertyToID("spiky_gradient_first_term");
    private readonly int _viscosityLaplacianFirstTermID = Shader.PropertyToID("viscosity_laplacian_first_term");
    
    private readonly int _stiffnessID = Shader.PropertyToID("stiffness");
    private readonly int _restDensityID = Shader.PropertyToID("rest_density");
    private readonly int _viscosityID = Shader.PropertyToID("viscosity");
    
    private readonly int _positionsBufferID = Shader.PropertyToID("positions");
    private readonly int _predictedPositionsBufferID = Shader.PropertyToID("predicted_positions");
    private readonly int _velocitiesBufferID = Shader.PropertyToID("velocities");
    private readonly int _densitiesBufferID = Shader.PropertyToID("densities");
    
    private readonly int _spatialKeysBufferID = Shader.PropertyToID("spatial_keys");
    private readonly int _spatialIndicesBufferID = Shader.PropertyToID("spatial_indices");
    private readonly int _spatialOffsetsBufferID = Shader.PropertyToID("spatial_offsets");

    private readonly int _auxiliaryPositionsBufferID = Shader.PropertyToID("auxiliary_positions");
    private readonly int _auxiliaryPredictedPositionsBufferID = Shader.PropertyToID("auxiliary_predicted_positions");
    private readonly int _auxiliaryVelocitiesBufferID = Shader.PropertyToID("auxiliary_velocities");
    
    #endregion

    #region GPU Related Fields
    
    private ComputeBuffer _positionsBuffer;
    private ComputeBuffer _predictedPositionsBuffer;
    private ComputeBuffer _velocitiesBuffer;
    private ComputeBuffer _densitiesBuffer;

    private ComputeBuffer _auxiliaryPositionsBuffer;
    private ComputeBuffer _auxiliaryPredictedPositionsBuffer;
    private ComputeBuffer _auxiliaryVelocitiesBuffer;

    private Texture2D _speedGradientTexture;

    private SpatialGrid _spatialGrid;

    private int _predictedAndKeysKernelID;
    private int _reorderKernelID;
    private int _copyAuxiliaryInMainKernelID;
    private int _densityKernelID;
    private int _deltaVelocityKernelID;

    #endregion

    [SerializeField] private uint[] preSortKeys;
    [SerializeField] private uint[] postSortKeys;
    [SerializeField] private uint[] offsets;
    [SerializeField] private uint[] indices;
    [SerializeField] private float[] densities;
    [SerializeField] private Vector3[] predictedPositions;
    [SerializeField] private Vector3[] positions;
    
    #region Initialization

    /// <summary>
    /// Initializes the entire simulation.
    /// </summary>
    private void Start()
    {
        // Linking the simulation compute shader
        _simulationComputeShader = Resources.Load<ComputeShader>("FluidSimulation");
        
        // Linking the fluid spawner
        _fluidSpawner = GetComponent<FluidSpawner>();
        
        // Initializing the particles
        InitializeParticles();
        
        // Initializing the variables related to optimization
        InitializeOptimization();
        
        // Initializing the kernel related variables
        UpdateKernelVariables();
        
        // Initializing the bounds of the simulation
        UpdateBoundsVariables();
        
        // Initializing the variables of the compute shader
        InitializeSimulationComputeShader();
        
        // Initializing the variables related to the rendering pipeline
        InitializeParticleShader();

        // Call to the debug function, used for debugging (duh)
        // preSortKeys = new uint[particlesAmount];
        // postSortKeys = new uint[particlesAmount];
        // offsets = new uint[particlesAmount];
        // indices = new uint[particlesAmount];
        // densities = new float[particlesAmount];
        // predictedPositions = new Vector3[particlesAmount];
        // positions = new Vector3[particlesAmount];
        // InternalDebug();
    }

    /// <summary>
    /// Wrapper function used to keep debugging stuff isolated
    /// </summary>
    private void InternalDebug()
    {
        // Dispatching the predicted position and spatial keys kernel
        _simulationComputeShader.Dispatch(_predictedAndKeysKernelID, Mathf.CeilToInt(_particlesAmount / 256.0f), 1, 1);
        
        // Debugging arrays
        _spatialGrid.SpatialKeysBuffer.GetData(preSortKeys);
        
        // Updating the spatial data
        _spatialGrid.UpdateSpatialLookup();
        
        // Debugging arrays
        _spatialGrid.SpatialKeysBuffer.GetData(postSortKeys);
        _spatialGrid.SpatialOffsetsBuffer.GetData(offsets);
        _spatialGrid.SpatialIndicesBuffer.GetData(indices);
        _predictedPositionsBuffer.GetData(predictedPositions);
        
        // Dispatching the reorder kernel
        _simulationComputeShader.Dispatch(_reorderKernelID, Mathf.CeilToInt(_particlesAmount / 256.0f), 1, 1);

        // Dispatching the copy auxiliary in main kernel
        _simulationComputeShader.Dispatch(_copyAuxiliaryInMainKernelID, 
            Mathf.CeilToInt(_particlesAmount / 256.0f), 1, 1);
        
        // Dispatching the density kernel
        _simulationComputeShader.Dispatch(_densityKernelID, _particlesAmount / 256 + 1, 1, 1);
        
        // Dispatching the delta velocities kernel
        _simulationComputeShader.Dispatch(_deltaVelocityKernelID, Mathf.CeilToInt(_particlesAmount / 256.0f), 1, 1);
        
        _densitiesBuffer.GetData(densities);
        _predictedPositionsBuffer.GetData(positions);
    }

    /// <summary>
    /// Initialized the variables related to the rendering of the simulation
    /// </summary>
    private void InitializeParticleShader()
    {
        // Assigning the buffer to the shader
        particleMaterial.SetBuffer(_positionsBufferID, _positionsBuffer);
        particleMaterial.SetBuffer(_velocitiesBufferID, _velocitiesBuffer);
        
        // Initializing the particles gradient texture
        _speedGradientTexture = new Texture2D(speedGradientTextureWidth, 1);
        
        // Setting the texture based on the gradient value
        for (var i = 0; i < speedGradientTextureWidth; i++)
        {
            // Mapping current value from [0, gradientTextureWidth] to [0, 1]
            var ratio = i / (float)(speedGradientTextureWidth - 1);
            
            // Obtaining the equivalent value in the gradient component setup
            var gradientColor = particleSpeedGradient.Evaluate(ratio);
            
            // Setting the pixel in the texture
            _speedGradientTexture.SetPixel(i, 0, gradientColor);
        }
        
        // Applying the changes to the texture
        _speedGradientTexture.Apply();
        
        // Passing the texture to the shader
        particleMaterial.SetTexture("_SpeedGradientTex", _speedGradientTexture);
    }
    
    /// <summary>
    /// Initializes the variables related to the kernels (smoothing functions).
    /// </summary>
    private void UpdateKernelVariables()
    {
        // Recomputing the squared kernel value
        _sqrKernelRadius = kernelRadius * kernelRadius;
        
        // Recomputing the volume of the kernels - 3D
        _poly6Normalization = 315.0f / (64.0f * Mathf.PI * Mathf.Pow(kernelRadius, 9));
        _spikyGradientNormalization = - 45.0f / (Mathf.PI * Mathf.Pow(kernelRadius, 6));
        _viscosityLaplacianNormalization = 45.0f / (Mathf.PI * Mathf.Pow(kernelRadius, 6));
    }

    /// <summary>
    /// Initializes the variables related to the bounds of the simulation
    /// </summary>
    private void UpdateBoundsVariables()
    {
        // Computing the maximum particles coordinates per axis
        _halfBounds = bounds / 2 - Vector3.one * particleRadius;
    }
    
    /// <summary>
    /// Initializes the variables related to the optimization.
    /// </summary>
    private void InitializeOptimization()
    {
        // CPU
        
        // Initializing the optimization related arrays
        _spatialLookup = new SpatialEntry[_particlesAmount];
        _lookupStartIndex = new int[_particlesAmount];
        
        // Creating the list of buckets
        _buckets = new List<SpatialEntry>[_particlesAmount];
        
        // Initializing the buckets
        for (var i = 0; i < _particlesAmount; i++)
        {
            _buckets[i] = new List<SpatialEntry>();
        }
        
        // GPU
        
        // Initializing the GPU implementation of the spatial grid
        _spatialGrid = new SpatialGrid(_particlesAmount);
    }

    /// <summary>
    /// Initializes the particles
    /// </summary>
    private void InitializeParticles()
    {
        // Initialzing the particles amount
        _particlesAmount = _fluidSpawner.GetParticlesAmount();
        
        // Initializing the particles related arrays
        _positions = new Vector3[_particlesAmount];
        _predictedPositions = new Vector3[_particlesAmount];
        _velocities = new Vector3[_particlesAmount];
        _densities = new float[_particlesAmount];
        
        // Initializing the gravity acceleration
        _gravityAcceleration = gravity * Vector3.down;
        
        // Setting the positions equal to the one provided by the fluid spawner
        _positions = _fluidSpawner.GetSpawnPositions();
    }

    /// <summary>
    /// Initializes the IDs of the ComputeShader's variables and initializes them.
    /// </summary>
    private void InitializeSimulationComputeShader()
    {
        // Initializing the compute buffer for the simulation on GPU
        _positionsBuffer = new ComputeBuffer(_particlesAmount, 3 * sizeof(float));
        _auxiliaryPositionsBuffer = new ComputeBuffer(_particlesAmount, 3 * sizeof(float));
        _predictedPositionsBuffer = new ComputeBuffer(_particlesAmount, 3 * sizeof(float));
        _auxiliaryPredictedPositionsBuffer = new ComputeBuffer(_particlesAmount, 3 * sizeof(float));
        _velocitiesBuffer = new ComputeBuffer(_particlesAmount, 3 * sizeof(float));
        _auxiliaryVelocitiesBuffer = new ComputeBuffer(_particlesAmount, 3 * sizeof(float));
        _densitiesBuffer = new ComputeBuffer(_particlesAmount, sizeof(float));
        
        // Filling the buffer
        _positionsBuffer.SetData(_positions);
        _predictedPositionsBuffer.SetData(_predictedPositions);
        _velocitiesBuffer.SetData(_velocities);
        _densitiesBuffer.SetData(_densities);
        
        // Setting the particles amount
        _simulationComputeShader.SetInt("particles_amount", _particlesAmount);
        _simulationComputeShader.SetFloat("particle_mass", particleMass);
        
        // Caching the reference ID of the ComputeShader's kernels
        _reorderKernelID = _simulationComputeShader.FindKernel("reorder_buffers");
        _copyAuxiliaryInMainKernelID = _simulationComputeShader.FindKernel("copy_auxiliary_in_main");
        _predictedAndKeysKernelID = _simulationComputeShader.FindKernel("calculate_predicted_and_keys");
        _densityKernelID = _simulationComputeShader.FindKernel("calculate_density");
        _deltaVelocityKernelID = _simulationComputeShader.FindKernel("calculate_delta_velocity");
        
        // Setting the buffers for the predicted positions kernel
        _simulationComputeShader.SetBuffer(_predictedAndKeysKernelID, _positionsBufferID, _positionsBuffer);
        _simulationComputeShader.SetBuffer(_predictedAndKeysKernelID, _predictedPositionsBufferID, _predictedPositionsBuffer);
        _simulationComputeShader.SetBuffer(_predictedAndKeysKernelID, _velocitiesBufferID, _velocitiesBuffer);
        _simulationComputeShader.SetBuffer(_predictedAndKeysKernelID, _spatialKeysBufferID, _spatialGrid.SpatialKeysBuffer);
        
        // Setting the buffers for the reorder kernel
        _simulationComputeShader.SetBuffer(_reorderKernelID, _spatialIndicesBufferID, _spatialGrid.SpatialIndicesBuffer);
        _simulationComputeShader.SetBuffer(_reorderKernelID, _positionsBufferID, _positionsBuffer);
        _simulationComputeShader.SetBuffer(_reorderKernelID, _predictedPositionsBufferID, _predictedPositionsBuffer);
        _simulationComputeShader.SetBuffer(_reorderKernelID, _velocitiesBufferID, _velocitiesBuffer);
        _simulationComputeShader.SetBuffer(_reorderKernelID, _auxiliaryPositionsBufferID, _auxiliaryPositionsBuffer);
        _simulationComputeShader.SetBuffer(_reorderKernelID, _auxiliaryPredictedPositionsBufferID, _auxiliaryPredictedPositionsBuffer);
        _simulationComputeShader.SetBuffer(_reorderKernelID, _auxiliaryVelocitiesBufferID, _auxiliaryVelocitiesBuffer);
        
        // Setting the buffers for the copy auxiliary in main kernel
        _simulationComputeShader.SetBuffer(_copyAuxiliaryInMainKernelID, _positionsBufferID, _positionsBuffer);
        _simulationComputeShader.SetBuffer(_copyAuxiliaryInMainKernelID, _predictedPositionsBufferID, _predictedPositionsBuffer);
        _simulationComputeShader.SetBuffer(_copyAuxiliaryInMainKernelID, _velocitiesBufferID, _velocitiesBuffer);
        _simulationComputeShader.SetBuffer(_copyAuxiliaryInMainKernelID, _auxiliaryPositionsBufferID, _auxiliaryPositionsBuffer);
        _simulationComputeShader.SetBuffer(_copyAuxiliaryInMainKernelID, _auxiliaryPredictedPositionsBufferID, _auxiliaryPredictedPositionsBuffer);
        _simulationComputeShader.SetBuffer(_copyAuxiliaryInMainKernelID, _auxiliaryVelocitiesBufferID, _auxiliaryVelocitiesBuffer);
        
        // Setting the buffers for the density kernel
        _simulationComputeShader.SetBuffer(_densityKernelID, _predictedPositionsBufferID, _predictedPositionsBuffer);
        _simulationComputeShader.SetBuffer(_densityKernelID, _densitiesBufferID, _densitiesBuffer);
        _simulationComputeShader.SetBuffer(_densityKernelID, _spatialKeysBufferID, _spatialGrid.SpatialKeysBuffer);
        _simulationComputeShader.SetBuffer(_densityKernelID, _spatialIndicesBufferID, _spatialGrid.SpatialIndicesBuffer);
        _simulationComputeShader.SetBuffer(_densityKernelID, _spatialOffsetsBufferID, _spatialGrid.SpatialOffsetsBuffer);
        
        // Setting the buffers for the delta velocity kernel
        _simulationComputeShader.SetBuffer(_deltaVelocityKernelID, _positionsBufferID, _positionsBuffer);
        _simulationComputeShader.SetBuffer(_deltaVelocityKernelID, _predictedPositionsBufferID, _predictedPositionsBuffer);
        _simulationComputeShader.SetBuffer(_deltaVelocityKernelID, _velocitiesBufferID, _velocitiesBuffer);
        _simulationComputeShader.SetBuffer(_deltaVelocityKernelID, _densitiesBufferID, _densitiesBuffer);
        _simulationComputeShader.SetBuffer(_deltaVelocityKernelID, _spatialKeysBufferID, _spatialGrid.SpatialKeysBuffer);
        _simulationComputeShader.SetBuffer(_deltaVelocityKernelID, _spatialIndicesBufferID, _spatialGrid.SpatialIndicesBuffer);
        _simulationComputeShader.SetBuffer(_deltaVelocityKernelID, _spatialOffsetsBufferID, _spatialGrid.SpatialOffsetsBuffer);
        
        UpdateComputeShaderVariables();
    }

    /// <summary>
    /// Updates the simulation compute shaders variables to the match the Unity's Editor.
    /// </summary>
    private void UpdateComputeShaderVariables()
    {
        // Setting ComputeShader's vectors
        _simulationComputeShader.SetVector(_halfBoundsID, _halfBounds);
        
        // Setting ComputeShader's matrices
        _simulationComputeShader.SetMatrix(_localToWorldMatrixID, transform.localToWorldMatrix);
        _simulationComputeShader.SetMatrix(_worldToLocalMatrixID, transform.worldToLocalMatrix);
        
        // Setting ComputeShader's floats
        _simulationComputeShader.SetFloat(_collisionDampingID, collisionDamping);
        
        _simulationComputeShader.SetFloat(_gravityID, gravity);
        
        _simulationComputeShader.SetFloat(_kernelRadiusID, kernelRadius);
        _simulationComputeShader.SetFloat(_sqrKernelRadiusID, _sqrKernelRadius);
        
        _simulationComputeShader.SetFloat(_poly6NormalizationID, _poly6Normalization);
        _simulationComputeShader.SetFloat(_spikyGradientFirstTermID, _spikyGradientNormalization);
        _simulationComputeShader.SetFloat(_viscosityLaplacianFirstTermID, _viscosityLaplacianNormalization);
        
        _simulationComputeShader.SetFloat(_restDensityID, restDensity);
        _simulationComputeShader.SetFloat(_stiffnessID, stiffness);
        _simulationComputeShader.SetFloat(_viscosityID, viscosity);
    }
    
    #endregion

    #region Unity Callback Functions

    /// <summary>
    /// Render the bounds of the simulation.
    /// </summary>
    private void OnDrawGizmos()
    {
        // If the simulation is running, return
        if (Application.isPlaying)
            return;
        
        // Drawing the simulation bounding box
        Gizmos.color = Color.white;
        Gizmos.DrawWireCube(new Vector3(0, 0, 0), bounds);
    }

    /// <summary>
    /// Clean-up of the ComputeBuffers used for rendering
    /// </summary>
    private void OnDestroy()
    {
        // Releasing Compute Buffers to prevent memory leaks
        _positionsBuffer?.Release();
        _predictedPositionsBuffer?.Release();
        _velocitiesBuffer?.Release();
        _densitiesBuffer?.Release();
        _auxiliaryPositionsBuffer?.Release();
        _auxiliaryPredictedPositionsBuffer?.Release();
        _auxiliaryVelocitiesBuffer?.Release();
        
        // Releasing the buffers allocated by the spatial grid
        _spatialGrid?.Release();
    }

    /// <summary>
    /// Updates the underlying simulation.
    /// </summary>
    private void Update()
    {
        // Updating the variables related to the kernels
        UpdateKernelVariables();
        
        // Updating the variables related to the bounds
        UpdateBoundsVariables();
        
        // Setting all the constants, in case they were changed via GUI
        UpdateComputeShaderVariables();
        
        // Computing a smaller delta time per substep
        var deltaTime = Mathf.Min(Time.deltaTime, 1 / 60.0f);
        
        // Calculating the delta time per simulation step
        var stepDeltaTime = deltaTime / simulationSubSteps;
        
        // Setting the delta time in the GPU
        _simulationComputeShader.SetFloat(_deltaTimeID, stepDeltaTime);
        
        // Running multiple simulation steps 
        for (var i = 0; i < simulationSubSteps; i++)
        {
            // Dispatching the predicted position kernel
            _simulationComputeShader.Dispatch(_predictedAndKeysKernelID, 
                Mathf.CeilToInt(_particlesAmount / 256.0f), 1, 1);
            
            // Updating the spatial data
            _spatialGrid.UpdateSpatialLookup();
            
            // Dispatching the reorder kernel
            _simulationComputeShader.Dispatch(_reorderKernelID, Mathf.CeilToInt(_particlesAmount / 256.0f), 1, 1);

            // Dispatching the copy auxiliary in main kernel
            _simulationComputeShader.Dispatch(_copyAuxiliaryInMainKernelID, 
                Mathf.CeilToInt(_particlesAmount / 256.0f), 1, 1);
            
            // Dispatching the density kernel
            _simulationComputeShader.Dispatch(_densityKernelID, Mathf.CeilToInt(_particlesAmount / 256.0f), 1, 1);
            
            // Dispatching the delta velocities kernel
            _simulationComputeShader.Dispatch(_deltaVelocityKernelID, Mathf.CeilToInt(_particlesAmount / 256.0f), 1, 1);
        }

        return;
        
        // Running multiple simulation steps 
        for (var i = 0; i < simulationSubSteps; i++)
        {
            EulerSimulationStep(Time.fixedDeltaTime / simulationSubSteps);
        }
    }

    /// <summary>
    /// Renders the simulation.
    /// </summary>
    private void LateUpdate()
    {
        // If the simulation is not to be rendered, return
        if (!renderSimulation)
            return;
        
        // Updating the uniform Radius on GPU
        particleMaterial.SetFloat(_particleRadiusID, particleRadius);

        // Updating position buffer on GPU
        // _positionsBuffer.SetData(_positions);
        // _velocitiesBuffer.SetData(_velocities);
        
        // Setting the rendering parameters
        var renderingParameters = new RenderParams(particleMaterial);
        
        // Render the particles
        Graphics.RenderPrimitives(renderingParameters, MeshTopology.Points, _particlesAmount);
    }
    
    #endregion

    #region Simulation
    
    /// <summary>
    /// Execute a simulation step with the given delta time using Euler based integration.
    /// </summary>
    /// <param name="deltaTime">The delta time from the previous simulation step.</param>
    private void EulerSimulationStep(float deltaTime)
    {
        // Computing the gravity vector
        _gravityAcceleration = gravity * deltaTime * Vector2.down;
        
        // Applying gravity and predict next positions
        Parallel.For(0, _particlesAmount, currentIndex => 
        {
            // Applying gravity to the particle velocity
            _velocities[currentIndex] += _gravityAcceleration;
            
            // Predicting the next particle position
            _predictedPositions[currentIndex] = _positions[currentIndex] + _velocities[currentIndex] * deltaTime;
        });
        
        // Updating the spatial lookup array
        UpdateSpatialLookup();
        
        // Calculating the density at particles position
        Parallel.For(0, _particlesAmount, currentIndex =>
        {
            // Computing the density at particle position
            CalculateDensity(currentIndex);
        });

        // Adding the current acceleration to the particle's velocity
        Parallel.For(0, _particlesAmount, currentIndex =>
        {
            // Adding the acceleration to the particle velocity
            _velocities[currentIndex] += CalculateAcceleration(currentIndex) * deltaTime;
        });
        
        // Updating the particles position and resolving collision
        Parallel.For(0, _particlesAmount, currentIndex =>
        {
            // Updating position
            _positions[currentIndex] += _velocities[currentIndex] * deltaTime;

            // Resolving collisions
            ResolveCollisions(ref _positions[currentIndex], ref _velocities[currentIndex]);
        });
    }

    #endregion

    #region Spatial Lookup

    /// <summary>
    /// Updates the spatial lookup data structure based on the current particles position.
    /// </summary>
    private void UpdateSpatialLookup()
    {
        // Creating the unordered spatial lookup array
        Parallel.For(0, _particlesAmount, currentIndex => {
            // Mapping the particle position to a cell in the grid
            var cellCoordinates = PositionToCellCoordinates(_predictedPositions[currentIndex]);
            
            // Calculating the cell key
            var cellKey = GetKeyFromHash(HashCell(cellCoordinates));
            
            // Locking the corresponding bucket to insert current particle
            lock (_buckets[cellKey])
            {
                _buckets[cellKey].Add(new SpatialEntry
                {
                    particleIndex = currentIndex,
                    cellKey = cellKey
                });
            }
            
            // Initializing the start index to infinity
            _lookupStartIndex[currentIndex] = int.MaxValue;
        });
        
        // Auxiliary variable used to flatten the buckets
        var spatialIndex = 0;
        
        // Flattening the buckets into a sorted array
        for (var i = 0; i < _particlesAmount; i++)
        {
            // Storing the current spatial index as the first with the current spatial key
            // Even if the bucket is empty it works because an empty bucket will never be accessed
            _lookupStartIndex[i] = spatialIndex;
            
            // Pushing all the values contained in the bucket into the array, resulting in a sorted array
            for (var j = 0; j < _buckets[i].Count; j++)
            {
                _spatialLookup[spatialIndex++] = _buckets[i][j];
            }
            
            // Clearing current bucket
            // TODO: Improve if more than 10000 particles
            _buckets[i].Clear();
        }
    }

    /// <summary>
    /// Given a position in a different coordinate system (ideally world space), returns the coordinates of the cell
    /// containing it.
    /// </summary>
    /// <param name="position">The 2D position to evaluate.</param>
    /// <returns>A tuple representing the coordinates of the cell containing the given position.</returns>
    private Vector3Int PositionToCellCoordinates(Vector3 position)
    {
        // Computing the coordinates with respect to the given radius
        var x = (int) (position.x / spatialGridCellSize);
        var y = (int) (position.y / spatialGridCellSize);
        var z = (int) (position.z / spatialGridCellSize);

        return new Vector3Int(x, y, z);
    }

    /// <summary>
    /// Given a 2D coordinate, returns their corresponding hash value.
    /// </summary>
    /// <param name="cellCoordinates">A tuple representing the cell coordinates.</param>
    /// <returns>The corresponding hash value of the coordinates.</returns>
    private int HashCell(Vector3Int cellCoordinates)
    {
        // Multiplying the coordinates by three arbitrary prime number
        var a = cellCoordinates.x * 73856093;
        var b = cellCoordinates.y * 19349663;
        var c = cellCoordinates.z * 83492791;
        
        // Returning the hashed value
        return a ^ b ^ c;
    }
    
    /// <summary>
    /// Given a hash value, returns its corresponding index in the spatial lookup array.
    /// </summary>
    /// <param name="hash">The hash value to map.</param>
    /// <returns>An index corresponding to the hash value.</returns>
    private int GetKeyFromHash(int hash)
    {
        return Mathf.Abs(hash) % _particlesAmount;
    }

    #endregion
    
    #region Collision Detection

    /// <summary>
    /// Given a reference to a particle's position and velocity, check collision with the simulation boundaries.
    /// </summary>
    /// <param name="position">A reference to the particle's position</param>
    /// <param name="velocity">A reference to the particle's velocity</param>
    private void ResolveCollisions(ref Vector3 position,ref Vector3 velocity)
    {
        // Case in which the particle collided with the bounds on X-axis
        if (Mathf.Abs(position.x) > _halfBounds.x)
        {
            // Setting the position x to the maximum
            position.x = _halfBounds.x * Mathf.Sign(position.x);
            // Flipping the velocity sign and applying a damping
            velocity.x *= -collisionDamping;
        }
        
        // Case in which the particle collided with the bounds on Y-axis
        if (Mathf.Abs(position.y) > _halfBounds.y)
        {
            // Setting the position y to the maximum
            position.y = _halfBounds.y * Mathf.Sign(position.y);
            // Flipping the velocity sign and applying a damping
            velocity.y *= -collisionDamping;
        }
        
        // Case in which the particle collided with the bound on Z-Axis
        if (Mathf.Abs(position.z) > _halfBounds.z)
        {
            // Setting the position y to the maximum
            position.z = _halfBounds.z * Mathf.Sign(position.z);
            // Flipping the velocity sign and applying a damping
            velocity.z *= -collisionDamping;
        }
    }

    #endregion
    
    #region Density
    
    /// <summary>
    /// Calculate the density at the position of the particle corresponding to given index.
    /// </summary>
    /// <param name="targetParticleIndex">The index of the particle to evaluate.</param>
    /// <returns>The density at the particle position</returns>
    private void CalculateDensity(int targetParticleIndex)
    {
        // Computing the coordinates of the grid cell containing the given position
        var cellCoordinates = PositionToCellCoordinates(_predictedPositions[targetParticleIndex]);
        
        // Initializing the density of the target particle
        var density = 0.0f;
        
        // Looping over all the cells belonging to the 3x3 square around the particle
        foreach (var offset in CellOffsets)
        {
            // Computing the key of the central cell + the current offset
            var currentCellKey = GetKeyFromHash(HashCell(new Vector3Int(
                cellCoordinates.x + offset.x, 
                cellCoordinates.y + offset.y, 
                cellCoordinates.z + offset.z)));
            
            // Extracting the index of the first element with current cell key
            var startIndex = _lookupStartIndex[currentCellKey];

            // Iterating all the particles belonging to the current cell
            for (var i = startIndex; i < _particlesAmount && _spatialLookup[i].cellKey == currentCellKey; i++)
            {
                // Extracting a reference to the current particle
                var currentParticleIndex = _spatialLookup[i].particleIndex;
                
                // Computing the distance between the current particle and the particle of interest
                var sqrDistance = (_predictedPositions[currentParticleIndex] 
                                   - _predictedPositions[targetParticleIndex]).sqrMagnitude;
            
                // Case in which the particle is outside the kernel radius or affecting itself
                if(sqrDistance >= _sqrKernelRadius)
                    continue;
            
                // Computing the influence of the current particle
                var influence = Poly6(sqrDistance);
            
                // Adding the result to the density
                density += particleMass * influence;
            }
        }
        
        _densities[targetParticleIndex] = density;
    }

    #endregion
    
    #region Pressure & Viscosity

    /// <summary>
    /// Computes the acceleration of the target particle, using the standard SPH formulation.
    /// </summary>
    /// <param name="targetParticleIndex"></param>
    /// <returns></returns>
    private Vector3 CalculateAcceleration(int targetParticleIndex)
    {
        // Computing the coordinates of the grid cell containing the given position
        var cellCoordinates = PositionToCellCoordinates(_predictedPositions[targetParticleIndex]);

        // Initializing the density
        var force = Vector3.zero;

        // Looping over all the cells belonging to the 3x3 square around the particle
        foreach (var offset in CellOffsets)
        {
            // Computing the key of the central cell + the current offset
            var currentCellKey = GetKeyFromHash(HashCell(new Vector3Int(
                cellCoordinates.x + offset.x, 
                cellCoordinates.y + offset.y, 
                cellCoordinates.z + offset.z)));

            // Extracting the index of the first element with current cell key
            var startIndex = _lookupStartIndex[currentCellKey];

            // Iterating all the particles belonging to the current cell
            for (var i = startIndex; i < _particlesAmount && _spatialLookup[i].cellKey == currentCellKey; i++)
            {
                // Extracting the index of the current particles
                var currentParticleIndex = _spatialLookup[i].particleIndex;
                
                // Skipping calculation for the particle by itself
                if (currentParticleIndex == targetParticleIndex)
                    continue;
                
                // Computing the vector from target particle to current particle
                var targetToCurrent =  _predictedPositions[targetParticleIndex] 
                                       - _predictedPositions[currentParticleIndex];
                
                // Case in which the particle is outside the kernel radius
                if(targetToCurrent.sqrMagnitude >= _sqrKernelRadius)
                    continue;
                
                // Computing the distance
                var distance = targetToCurrent.magnitude;
                
                // Computing the influence of the current particle
                var viscosityInfluence = ViscosityLaplacian(distance);
                
                // Computing the difference between velocities
                var velocitiesDifference = _velocities[currentParticleIndex] - _velocities[targetParticleIndex];
                
                // Adding the result to the force
                force += particleMass * viscosityInfluence * viscosity * velocitiesDifference 
                         / _densities[currentParticleIndex];
                
                // Computing the direction between the current particle and the point of interest
                var pressureDirection = distance == 0 ?
                    new Vector3(1, 0, 0)
                    : targetToCurrent / distance;
            
                // Evaluating the slope of the kernel at the given distance
                var kernelSlope = SpikyGradient(distance);
            
                // Computing the shared pressure force
                var sharedPressure = CalculateSharedPressure(_densities[currentParticleIndex],
                    _densities[targetParticleIndex]);
                
                // Computing the applied force using both pressure and viscosity
                force -= particleMass * sharedPressure * kernelSlope * pressureDirection 
                         / _densities[currentParticleIndex];
            }
        }
        
        return force / _densities[targetParticleIndex] + _gravityAcceleration;
    }
    
    /// <summary>
    /// Convert the provided density value into pressure using SPH formula
    /// </summary>
    /// <param name="density">The provided density</param>
    /// <returns></returns>
    private float ConvertDensityToPressure(float density)
    {
        // Computing the pressure force using Tait equation
        return stiffness * (Mathf.Pow(density / restDensity, 4) - 1);
    }

    /// <summary>
    /// Auxiliary function used to compute the pressure between two particles given the density at their position.
    /// </summary>
    /// <param name="firstDensity">Density at the first particle's position.</param>
    /// <param name="secondDensity">Density at the second particle's position.</param>
    /// <returns>The pressure force between the particles.</returns>
    private float CalculateSharedPressure(float firstDensity, float secondDensity)
    {
        // Converting both densities to pressure 
        var firstPressure = ConvertDensityToPressure(firstDensity);
        var secondPressure = ConvertDensityToPressure(secondDensity);

        // Computing the pressure that respects Newton's Third Law of Motion (as explained by Sebastian Lague)
        return (firstPressure + secondPressure) / 2;
    }

    #endregion

    #region Kernels
    
    /// <summary>
    /// Computes the gradient of the Spiky kernel at the given distance
    /// </summary>
    /// <param name="distance">The distance</param>
    /// <returns>The gradient of the kernel at given distance.</returns>
    private float SpikyGradient(float distance)
    {
        // Computing the gradient 
        return _spikyGradientNormalization * Mathf.Pow(distance - kernelRadius, 2);
    }

    /// <summary>
    /// Computes the influence at the given squared distance using the kernel: (radius^2 - distance^2)^3
    /// </summary>
    /// <param name="sqrDistance">The squared distance</param>
    /// <returns>The influence at the given distance</returns>
    private float Poly6(float sqrDistance)
    {
        // Computing the influence at given distance
        return _poly6Normalization * Mathf.Pow(_sqrKernelRadius - sqrDistance, 3);
    }
    
    /// <summary>
    /// Computes the laplacian of the Viscosity kernel at the given distance
    /// </summary>
    /// <param name="distance">The distance to evaluate.</param>
    /// <returns>The laplacian of the kernel at given distance.</returns>
    private float ViscosityLaplacian(float distance)
    {
        //Computing the laplacian
        return _viscosityLaplacianNormalization * (kernelRadius - distance);
    }

    #endregion
}