using System.Collections.Generic;
using System.Threading.Tasks;
using UnityEngine;
using UnityEngine.Serialization;
using UnityEngine.XR;
using Random = UnityEngine.Random;
using Vector2 = UnityEngine.Vector2;
using Vector3 = UnityEngine.Vector3;

public class FluidSimulation : MonoBehaviour
{
    [Header("Particle")]
    public float particleRadius;
    public float particleMass;
    [Range(0.0f, 1.0f)] public float collisionDamping;
    
    [Header("Kernel")]
    public float kernelRadius;
    
    [Header("Bounds")]
    public Vector3 boundsSize;
    
    [Header("Simulation")]
    public float gravityIntensity;
    public Vector3 gravityDirection;
    public float externalForcesIntensity;
    public float restDensity;
    public float stiffness;
    public float viscosity;
    public float nearDensity;
    public uint simulationSubSteps;
    public float slowMotionCoefficient;
    
    [Header("Rendering")]
    public Material particleMaterial;
    public Gradient particleSpeedGradient;
    public int speedGradientTextureWidth;
    public bool renderSimulation;
    public RayMarchedFluid rayMarchedFluid;
    public RayMarchedGas rayMarchedGas;
    public RayMarchedNormals rayMarchedNormals;
    public DensitySlice densitySlice;
    
    [Header("Density Map")] 
    public Vector3Int densityMapSize;
    public bool updateDensityMap;

    [Header("Floor")] 
    public Transform floorTransform;
    
    #region CPU Related Fields
    
    // Reference to the simulation compute shader
    private ComputeShader _simulationComputeShader;
    
    // Reference to the particles spawner
    private FluidSpawner _fluidSpawner;
    // Density Map
    private RenderTexture _densityMap;
    
    // Fields related to user input
    private bool _isSimulationPlaying = true;
    private bool _isSimulationSlowMotion = false;
    
    // Particles related fields
    private int _particlesAmount;
    private Vector3[] _positions;
    private Vector3[] _predictedPositions;
    private Vector3[] _velocities;
    private float[] _densities;
    private float[] _nearDensities;
    private Vector3 _externalForcesDirection;
    
    // Bounds related fields
    private Vector3 _halfBoundsSize;
    
    // Kernel related fields
    private float _sqrKernelRadius;
    private float _poly6Normalization;
    private float _spikyGradientNormalization;
    private float _viscosityLaplacianNormalization;
    
    // Shaders related fields
    private readonly int _halfBoundsSizeID = Shader.PropertyToID("half_bounds_size");
    private readonly int _boundsSizeID = Shader.PropertyToID("bounds_size");
    private readonly int _collisionDampingID = Shader.PropertyToID("collision_damping");

    private readonly int _worldToLocalMatrixID = Shader.PropertyToID("world_to_local");
    private readonly int _localToWorldMatrixID = Shader.PropertyToID("local_to_world");

    private readonly int _deltaTimeID = Shader.PropertyToID("delta_time");
    
    private readonly int _densityMapID = Shader.PropertyToID("density_map");
    private readonly int _densityMapPropertyID = Shader.PropertyToID("_DensityMap");
    private readonly int _densityMapSizeID = Shader.PropertyToID("density_map_size");
    
    private readonly int _externalForcesAccelerationID = Shader.PropertyToID("external_forces_acceleration");
    
    private readonly int _particleRadiusID = Shader.PropertyToID("particle_radius");
    private readonly int _kernelRadiusID = Shader.PropertyToID("kernel_radius");
    private readonly int _sqrKernelRadiusID = Shader.PropertyToID("sqr_kernel_radius");
    
    private readonly int _poly6NormalizationID = Shader.PropertyToID("poly6_normalization");
    private readonly int _spikyGradientFirstTermID = Shader.PropertyToID("spiky_gradient_first_term");
    private readonly int _viscosityLaplacianFirstTermID = Shader.PropertyToID("viscosity_laplacian_first_term");

    private readonly int _particleMassID = Shader.PropertyToID("particle_mass");
    private readonly int _stiffnessID = Shader.PropertyToID("stiffness");
    private readonly int _restDensityID = Shader.PropertyToID("rest_density");
    private readonly int _viscosityID = Shader.PropertyToID("viscosity");
    private readonly int _nearDensityID = Shader.PropertyToID("near_density");
    
    private readonly int _positionsBufferID = Shader.PropertyToID("positions");
    private readonly int _predictedPositionsBufferID = Shader.PropertyToID("predicted_positions");
    private readonly int _velocitiesBufferID = Shader.PropertyToID("velocities");
    private readonly int _densitiesBufferID = Shader.PropertyToID("densities");
    private readonly int _nearDensitiesBufferID = Shader.PropertyToID("near_densities");
    
    private readonly int _spatialKeysBufferID = Shader.PropertyToID("spatial_keys");
    private readonly int _spatialIndicesBufferID = Shader.PropertyToID("spatial_indices");
    private readonly int _spatialOffsetsBufferID = Shader.PropertyToID("spatial_offsets");

    private readonly int _auxiliaryPositionsBufferID = Shader.PropertyToID("auxiliary_positions");
    private readonly int _auxiliaryPredictedPositionsBufferID = Shader.PropertyToID("auxiliary_predicted_positions");
    private readonly int _auxiliaryVelocitiesBufferID = Shader.PropertyToID("auxiliary_velocities");

    private readonly int _floorObjectToWorldID = Shader.PropertyToID("floor_object_to_world");
    private readonly int _floorWorldToObjectID = Shader.PropertyToID("floor_world_to_object");
    #endregion

    #region GPU Related Fields
    
    // Simulation Compute Shader's buffers
    private ComputeBuffer _positionsBuffer;
    private ComputeBuffer _predictedPositionsBuffer;
    private ComputeBuffer _velocitiesBuffer;
    private ComputeBuffer _densitiesBuffer;
    private ComputeBuffer _nearDensitiesBuffer;
    private ComputeBuffer _auxiliaryPositionsBuffer;
    private ComputeBuffer _auxiliaryPredictedPositionsBuffer;
    private ComputeBuffer _auxiliaryVelocitiesBuffer;
    
    // Particles material gradient texture
    private Texture2D _speedGradientTexture;
    
    // Reference to the spatial grid algorithm
    private SpatialGrid _spatialGrid;

    // Kernel IDs
    private int _predictedPositionKernelID;
    private int _spatialKeysKernelID;
    private int _reorderKernelID;
    private int _copyAuxiliaryInMainKernelID;
    private int _densitiesKernelID;
    private int _positionAndVelocityKernelID;
    private int _updateDensityMapKernelID;

    #endregion
    
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
        
        // Initializing the density map used in multiple rendering pipelines
        InitializeDensityMap();
        
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
    }

    /// <summary>
    /// Wrapper function used to keep debugging stuff isolated
    /// </summary>
    private void InternalDebug()
    {
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
        
        // Enable GPU instancing using the material
        particleMaterial.enableInstancing = true;
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
        _halfBoundsSize = boundsSize / 2;
    }
    
    /// <summary>
    /// Initializes the variables related to the optimization.
    /// </summary>
    private void InitializeOptimization()
    {
        // Initializing the GPU implementation of the spatial grid
        _spatialGrid = new SpatialGrid(_particlesAmount);
    }

    /// <summary>
    /// Initializes the particles
    /// </summary>
    private void InitializeParticles()
    {
        // Initializing the particles amount
        _particlesAmount = _fluidSpawner.GetParticlesAmount();
        
        // Initializing the particles related arrays
        _positions = new Vector3[_particlesAmount];
        _predictedPositions = new Vector3[_particlesAmount];
        _velocities = new Vector3[_particlesAmount];
        _densities = new float[_particlesAmount];
        _nearDensities = new float[_particlesAmount];
        
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
        _nearDensitiesBuffer = new ComputeBuffer(_particlesAmount, sizeof(float));
        
        // Filling the buffer
        _positionsBuffer.SetData(_positions);
        _predictedPositionsBuffer.SetData(_predictedPositions);
        _velocitiesBuffer.SetData(_velocities);
        _densitiesBuffer.SetData(_densities);
        
        // Setting not dynamic variables
        _simulationComputeShader.SetInt("particles_amount", _particlesAmount);
        _simulationComputeShader.SetFloat("particle_mass", particleMass);
        
        // Caching the reference ID of the ComputeShader's kernels
        _reorderKernelID = _simulationComputeShader.FindKernel("reorder_buffers");
        _copyAuxiliaryInMainKernelID = _simulationComputeShader.FindKernel("copy_auxiliary_in_main");
        _predictedPositionKernelID = _simulationComputeShader.FindKernel("update_predicted_positions");
        _spatialKeysKernelID = _simulationComputeShader.FindKernel("update_spatial_keys");
        _densitiesKernelID = _simulationComputeShader.FindKernel("update_densities");
        _positionAndVelocityKernelID = _simulationComputeShader.FindKernel("update_position_and_velocity");
        _updateDensityMapKernelID = _simulationComputeShader.FindKernel("update_density_map");
            
        // Setting the buffers for the predicted positions kernel
        _simulationComputeShader.SetBuffer(_predictedPositionKernelID, _positionsBufferID, _positionsBuffer);
        _simulationComputeShader.SetBuffer(_predictedPositionKernelID, _predictedPositionsBufferID, _predictedPositionsBuffer);
        _simulationComputeShader.SetBuffer(_predictedPositionKernelID, _velocitiesBufferID, _velocitiesBuffer);
        
        // Setting the buffers for the spatial keys kernel
        _simulationComputeShader.SetBuffer(_spatialKeysKernelID, _predictedPositionsBufferID, _predictedPositionsBuffer);
        _simulationComputeShader.SetBuffer(_spatialKeysKernelID, _spatialKeysBufferID, _spatialGrid.SpatialKeysBuffer);
        
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
        _simulationComputeShader.SetBuffer(_densitiesKernelID, _predictedPositionsBufferID, _predictedPositionsBuffer);
        _simulationComputeShader.SetBuffer(_densitiesKernelID, _densitiesBufferID, _densitiesBuffer);
        _simulationComputeShader.SetBuffer(_densitiesKernelID, _nearDensitiesBufferID, _nearDensitiesBuffer);
        _simulationComputeShader.SetBuffer(_densitiesKernelID, _spatialKeysBufferID, _spatialGrid.SpatialKeysBuffer);
        _simulationComputeShader.SetBuffer(_densitiesKernelID, _spatialIndicesBufferID, _spatialGrid.SpatialIndicesBuffer);
        _simulationComputeShader.SetBuffer(_densitiesKernelID, _spatialOffsetsBufferID, _spatialGrid.SpatialOffsetsBuffer);
        
        // Setting the buffers for the delta velocity kernel
        _simulationComputeShader.SetBuffer(_positionAndVelocityKernelID, _positionsBufferID, _positionsBuffer);
        _simulationComputeShader.SetBuffer(_positionAndVelocityKernelID, _predictedPositionsBufferID, _predictedPositionsBuffer);
        _simulationComputeShader.SetBuffer(_positionAndVelocityKernelID, _velocitiesBufferID, _velocitiesBuffer);
        _simulationComputeShader.SetBuffer(_positionAndVelocityKernelID, _densitiesBufferID, _densitiesBuffer);
        _simulationComputeShader.SetBuffer(_positionAndVelocityKernelID, _nearDensitiesBufferID, _nearDensitiesBuffer);
        _simulationComputeShader.SetBuffer(_positionAndVelocityKernelID, _spatialKeysBufferID, _spatialGrid.SpatialKeysBuffer);
        _simulationComputeShader.SetBuffer(_positionAndVelocityKernelID, _spatialIndicesBufferID, _spatialGrid.SpatialIndicesBuffer);
        _simulationComputeShader.SetBuffer(_positionAndVelocityKernelID, _spatialOffsetsBufferID, _spatialGrid.SpatialOffsetsBuffer);
        
        // Setting the texture for the density map kernel
        _simulationComputeShader.SetTexture(_updateDensityMapKernelID, _densityMapID, _densityMap);
        _simulationComputeShader.SetBuffer(_updateDensityMapKernelID, _predictedPositionsBufferID, _predictedPositionsBuffer);
        _simulationComputeShader.SetBuffer(_updateDensityMapKernelID, _spatialKeysBufferID, _spatialGrid.SpatialKeysBuffer);
        _simulationComputeShader.SetBuffer(_updateDensityMapKernelID, _spatialIndicesBufferID, _spatialGrid.SpatialIndicesBuffer);
        _simulationComputeShader.SetBuffer(_updateDensityMapKernelID, _spatialOffsetsBufferID, _spatialGrid.SpatialOffsetsBuffer);
        
        UpdateComputeShaderVariables();
    }

    /// <summary>
    /// Updates the simulation compute shaders variables to the match the Unity's Editor.
    /// </summary>
    private void UpdateComputeShaderVariables()
    {
        // Setting ComputeShader's vectors
        _simulationComputeShader.SetVector(_halfBoundsSizeID, _halfBoundsSize);
        _simulationComputeShader.SetVector(_boundsSizeID, boundsSize);
        _simulationComputeShader.SetInts(_densityMapSizeID, densityMapSize.x, densityMapSize.y, densityMapSize.z);
        
        // Setting ComputeShader's matrices
        _simulationComputeShader.SetMatrix(_localToWorldMatrixID, transform.localToWorldMatrix);
        _simulationComputeShader.SetMatrix(_worldToLocalMatrixID, transform.worldToLocalMatrix);
        
        // Setting ComputeShader's floats
        _simulationComputeShader.SetFloat(_collisionDampingID, collisionDamping);
        _simulationComputeShader.SetFloat(_particleMassID, particleMass);
        
        _simulationComputeShader.SetFloat(_kernelRadiusID, kernelRadius);
        _simulationComputeShader.SetFloat(_sqrKernelRadiusID, _sqrKernelRadius);
        
        _simulationComputeShader.SetFloat(_poly6NormalizationID, _poly6Normalization);
        _simulationComputeShader.SetFloat(_spikyGradientFirstTermID, _spikyGradientNormalization);
        _simulationComputeShader.SetFloat(_viscosityLaplacianFirstTermID, _viscosityLaplacianNormalization);
        
        _simulationComputeShader.SetFloat(_restDensityID, restDensity);
        _simulationComputeShader.SetFloat(_stiffnessID, stiffness);
        _simulationComputeShader.SetFloat(_viscosityID, viscosity);
        _simulationComputeShader.SetFloat(_nearDensityID, nearDensity);
    }

    /// <summary>
    /// Initializes the density map used by multiple rendering materials.
    /// </summary>
    private void InitializeDensityMap()
    {
        // Initializing the density map
        _densityMap = new RenderTexture(densityMapSize.x, densityMapSize.y, 0, RenderTextureFormat.RFloat)
        {
            // Signaling to Unity that this is a 3D texture
            dimension = UnityEngine.Rendering.TextureDimension.Tex3D,
            // Setting the texture depth 
            volumeDepth = densityMapSize.z,
            // Allowing random write on the texture (done by the Compute Shader)
            enableRandomWrite = true,
            // Setting the texture filter
            filterMode = FilterMode.Bilinear,
            // Setting the wrap mode
            wrapMode = TextureWrapMode.Clamp
        };

        // Creating the density map
        _densityMap.Create();

        // Assigning the map to every material that uses it
        densitySlice.material.SetTexture(_densityMapPropertyID, _densityMap);
        rayMarchedGas.material.SetTexture(_densityMapPropertyID, _densityMap);
        rayMarchedFluid.material.SetTexture(_densityMapPropertyID, _densityMap);
        rayMarchedNormals.material.SetTexture(_densityMapPropertyID, _densityMap);
    }

    /// <summary>
    /// Updates the global variables shared between multiple shaders.
    /// </summary>
    private void UpdateShaderGlobalVariables()
    {
        // Updating the floor global matrix
        Shader.SetGlobalMatrix(_floorWorldToObjectID, floorTransform.worldToLocalMatrix);
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
        Gizmos.DrawWireCube(new Vector3(0, 0, 0), boundsSize);
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
        _nearDensitiesBuffer?.Release();
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
        // Handling possible user's inputs
        HandleInputs();
        
        // If the simulation is not running, do no execute simulation step
        if (!_isSimulationPlaying)
            return;
        
        // Executing a simulation step
        ExecuteSimulationStep();
    }

    /// <summary>
    /// Executes a simulation step
    /// </summary>
    private void ExecuteSimulationStep()
    { 
        // Updating the scale of the renderings
        rayMarchedFluid.transform.localScale = boundsSize;
        rayMarchedGas.transform.localScale = boundsSize;
        rayMarchedNormals.transform.localScale = boundsSize;
        
        // Updating the variables related to the kernels
        UpdateKernelVariables();
        
        // Updating the variables related to the bounds
        UpdateBoundsVariables();
        
        // Setting all the constants, in case they were changed via GUI
        UpdateComputeShaderVariables();
        
        // Updating the shaders global variables
        UpdateShaderGlobalVariables();
        
        // Computing a smaller delta time per substep
        var deltaTime = Mathf.Min(Time.deltaTime, 1 / 60.0f);
        
        // Calculating the delta time per simulation step
        var stepDeltaTime = deltaTime / simulationSubSteps;
        
        // If the simulation is in slow motion, multiply the delta time by the slow motion coefficient
        if (_isSimulationSlowMotion)
            stepDeltaTime *= slowMotionCoefficient;
        
        // Setting the delta time in the GPU
        _simulationComputeShader.SetFloat(_deltaTimeID, stepDeltaTime);
        _simulationComputeShader.SetVector(_externalForcesAccelerationID, 
            stepDeltaTime * gravityIntensity * gravityDirection 
            + stepDeltaTime * externalForcesIntensity * _externalForcesDirection);
        
        // Declaring the const amount of threads for linear scan
        const float linearThreadsAmount = 512.0f;
        
        // Running multiple simulation steps 
        for (var i = 0; i < simulationSubSteps; i++)
        {
            // Dispatching the predicted position kernel
            _simulationComputeShader.Dispatch(_predictedPositionKernelID, 
                Mathf.CeilToInt(_particlesAmount / linearThreadsAmount), 1, 1);
            
            // Dispatching the spatial keys kernel
            _simulationComputeShader.Dispatch(_spatialKeysKernelID, 
                Mathf.CeilToInt(_particlesAmount / linearThreadsAmount), 1, 1);
            
            // Updating the spatial data
            _spatialGrid.Update();
    
            // Dispatching the reorder kernel
            _simulationComputeShader.Dispatch(_reorderKernelID, 
                Mathf.CeilToInt(_particlesAmount / linearThreadsAmount), 1, 1);

            // Dispatching the copy auxiliary in main kernel
            _simulationComputeShader.Dispatch(_copyAuxiliaryInMainKernelID, 
                Mathf.CeilToInt(_particlesAmount / linearThreadsAmount), 1, 1);
            
            // Dispatching the density kernel
            _simulationComputeShader.Dispatch(_densitiesKernelID, 
                Mathf.CeilToInt(_particlesAmount / linearThreadsAmount), 1, 1);
            
            // Dispatching the delta velocities kernel
            _simulationComputeShader.Dispatch(_positionAndVelocityKernelID,
                Mathf.CeilToInt(_particlesAmount / linearThreadsAmount), 1, 1);
        }

        // Updating the density map
        if(updateDensityMap)
            _simulationComputeShader.Dispatch(_updateDensityMapKernelID,
                Mathf.CeilToInt(densityMapSize.x / 8.0f), 
                Mathf.CeilToInt(densityMapSize.y / 8.0f), 
                Mathf.CeilToInt(densityMapSize.z / 8.0f));
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

        // Setting the rendering parameters
        var renderingParameters = new RenderParams(particleMaterial);
        
        // Render the particles
        Graphics.RenderPrimitives(renderingParameters, MeshTopology.Points, _particlesAmount);
    }
    
    #endregion
    
    # region Input

    /// <summary>
    /// Handles the input pressed by the user.
    /// </summary>
    private void HandleInputs()
    {
        // Pausing or starting the simulation
        if (Input.GetKeyDown(KeyCode.Space))
        {
            _isSimulationPlaying = !_isSimulationPlaying;
        }
        
        // Resetting the simulation
        if (Input.GetKeyDown(KeyCode.R))
        {
            // Resetting the particles data
            _positionsBuffer.SetData(_positions);
            _predictedPositionsBuffer.SetData(_predictedPositions);
            _densitiesBuffer.SetData(_densities);
            _nearDensitiesBuffer.SetData(_nearDensities);
            _velocitiesBuffer.SetData(_velocities);
            
            // Updating the data
            ExecuteSimulationStep();
        }
        
        // Disabling or enabling slow motion
        if (Input.GetKeyDown(KeyCode.F))
        {
            _isSimulationSlowMotion = !_isSimulationSlowMotion;
        }
        
        // Resetting the external forces direction
        _externalForcesDirection = Vector3.zero;
        
        // Adding the contribution of every possible external forces direction
        if (Input.GetKey(KeyCode.UpArrow))
            _externalForcesDirection += Vector3.forward;
        
        if (Input.GetKey(KeyCode.DownArrow))
            _externalForcesDirection += Vector3.back;
        
        if (Input.GetKey(KeyCode.RightArrow))
            _externalForcesDirection += Vector3.right;
        
        if (Input.GetKey(KeyCode.LeftArrow))
            _externalForcesDirection += Vector3.left;
        
        // Normalizing the direction
        _externalForcesDirection.Normalize();
    }
    
    #endregion
}