using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using UnityEngine;
using UnityEngine.Serialization;
using Vector2 = UnityEngine.Vector2;
using Vector3 = UnityEngine.Vector3;

public class FluidSimulation : MonoBehaviour
{
    [Header("Particle")]
    [SerializeField] public float particleRadius;
    [SerializeField] public float particleMass;
    [SerializeField] public float collisionDamping;

    [Header("Particles")]
    [SerializeField] public int particlesAmount;
    [SerializeField] public float particlesSpacing;
    [SerializeField] public ParticlesInitialization initializationAlgorithm;
    
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
    [SerializeField] public IntegrationMethod integrationMethod;
    [SerializeField] public ComputeShader simulationComputeShader;

    [Header("Optimization")] 
    [SerializeField] public float spatialGridCellSize;
    
    [Header("Rendering")]
    [SerializeField] public Material particleMaterial;
    [SerializeField] public Gradient particleSpeedGradient;
    [SerializeField] public int speedGradientTextureWidth;
    [SerializeField] public bool renderSimulation;
    

    #region CPU Related Fields

    [Serializable]
    private struct Particle
    {
        public Vector3 position;
        public Vector3 predictedPosition;
        public Vector3 velocity;
        public Vector3 acceleration;
        public float density;
    }
    
    public enum ParticlesInitialization
    {
        Grid,
        Random
    }

    public enum IntegrationMethod
    {
        Euler,
        LeapFrog
    }
    
    private struct SpatialEntry
    {
        public int particleIndex;
        public int cellKey;
    }
    
    // Particles related fields
    private Particle[] _particles;
    
    // Bounds related fields
    private Vector3 _halfBounds;
    
    // Simulation related fields
    private Vector3 _gravityAcceleration;
    
    // Data structure related fields
    private int _isPingActive;
    
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
    private int _halfBoundsID;
    private int _collisionDampingID;
    
    private int _gravityID;
    
    private int _poly6NormalizationID;
    private int _spikyGradientFirstTermID;
    private int _viscosityLaplacianFirstTermID;
    private int _kernelRadiusID;
    private int _sqrKernelRadiusID;

    private int _stiffnessID;
    private int _restDensityID;
    private int _viscosityID;

    private int _isPingActiveID;
    
    #endregion

    #region GPU Related Fields
    
    private ComputeBuffer _particlesBufferPing;
    private ComputeBuffer _particlesBufferPong;

    private Texture2D _speedGradientTexture;

    private int _predictedPositionKernelID;
    private int _densityKernelID;
    private int _deltaVelocityKernelID;

    #endregion

    #region Initialization

    /// <summary>
    /// Initializes the entire simulation.
    /// </summary>
    private void Start()
    {
        // Initializing the variables related to optimization
        InitializeOptimization();
        
        // Initializing the kernel related variables
        UpdateKernelVariables();
        
        // Initializing the bounds of the simulation
        UpdateBoundsVariables();
        
        // Initializing the particles
        InitializeParticles();
        
        // Initializing the variables of the compute shader
        InitializeComputeShader();
        
        // Initializing the variables related to the rendering pipeline
        InitializeParticleShader();
    }

    /// <summary>
    /// Initialized the variables related to the rendering of the simulation
    /// </summary>
    private void InitializeParticleShader()
    {
        // Assigning the buffer to the shader
        particleMaterial.SetBuffer("particles_buffer_ping", _particlesBufferPing);
        particleMaterial.SetBuffer("particles_buffer_pong", _particlesBufferPing);
        
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
        
        // Recomputing the volume of the kernels - 2D
        // TODO: 3D conversion of the simulation
        // _poly6Normalization = 4.0f / (Mathf.PI * Mathf.Pow(kernelRadius, 8));
        // _spikyGradientNormalization = -30.0f / (Mathf.PI * Mathf.Pow(kernelRadius, 5));
        // _viscosityLaplacianNormalization = 20.0f / (Mathf.PI * Mathf.Pow(kernelRadius, 5));
        
        // Recomputing the volume of the kernels - 3D
        _poly6Normalization = 315.0f / (64.0f * Mathf.PI * Mathf.Pow(kernelRadius, 9));
        _spikyGradientNormalization = -45.0f / (Mathf.PI * Mathf.Pow(kernelRadius, 6));
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
        // Initializing the optimization related arrays
        _spatialLookup = new SpatialEntry[particlesAmount];
        _lookupStartIndex = new int[particlesAmount];
        
        // Creating the list of buckets
        _buckets = new List<SpatialEntry>[particlesAmount];
        
        // Initializing the buckets
        for (var i = 0; i < particlesAmount; i++)
        {
            _buckets[i] = new List<SpatialEntry>();
        }
    }

    /// <summary>
    /// Initializes the particles
    /// </summary>
    private void InitializeParticles()
    {
        // Initializing the particles related arrays
        _particles = new Particle[particlesAmount];
        
        // Initializing the gravity acceleration
        _gravityAcceleration = gravity * Vector3.down;
        
        // Initializing the compute buffer for the GPU
        _particlesBufferPing = new ComputeBuffer(particlesAmount, (3 * 4 + 1) * sizeof(float));
        _particlesBufferPong = new ComputeBuffer(particlesAmount, (3 * 4 + 1) * sizeof(float));
        
        // Initializing the particles position based on the algorithm selected via editor
        switch (initializationAlgorithm)
        {
            case ParticlesInitialization.Random:
                for (var i = 0; i < particlesAmount; i++)
                {
                    // Computing the new coordinates
                    var randomX = (UnityEngine.Random.value - 0.5f) * bounds.x;
                    var randomY = (UnityEngine.Random.value - 0.5f) * bounds.y;
                    var randomZ = (UnityEngine.Random.value - 0.5f) * bounds.z;
                    
                    // Assigning the new position
                    _particles[i].position = new Vector3(randomX, randomY, randomZ);
                }
                break;
            
            case ParticlesInitialization.Grid:
            default:
                // Initializing the auxiliary variables used to create the particles grid
                var particlesPerRow = (int)Mathf.Sqrt(particlesAmount);
                var particlesPerCol = (particlesAmount - 1) / particlesPerRow + 1;
                var spacing = particleRadius * 2 + particlesSpacing;
                
                // Initializing the position of each particle
                for (var i = 0; i < particlesAmount; i++)
                {
                    // Computing the indexes of the current particle
                    var zIndex = i / (particlesPerRow * particlesPerRow);
                    var yIndex = (i / particlesPerRow) % particlesPerRow;
                    var xIndex = i % particlesPerRow;

                    // Computing the new coordinates
                    var x = (xIndex - particlesPerRow / 2.0f + 0.5f) * spacing;
                    var y = (yIndex - particlesPerRow / 2.0f + 0.5f) * spacing;
                    var z = (zIndex - particlesPerRow / 2.0f + 0.5f) * spacing;

                    _particles[i].position = new Vector3(x, y, z);
                }
                break;
        }
        
        // Initializing the particles acceleration if LeapFrog is enabled
        if (integrationMethod == IntegrationMethod.LeapFrog)
        {
            // Updating the spatial grid
            UpdateSpatialLookup();
            
            // Calculating the density at particles position
            Parallel.For(0, particlesAmount, currentIndex =>
            {
                // Computing the density at particle position
                CalculateDensity(ref _particles[currentIndex]);
            });
            
            // Calculating the acceleration pf particles
            Parallel.For(0, particlesAmount, currentIndex =>
            {
                _particles[currentIndex].acceleration = CalculateAcceleration(currentIndex);
            });
        }
        
        // Filling the buffer
        // TODO: 3D Conversion
        _particlesBufferPing.SetData(_particles);
        _particlesBufferPong.SetData(_particles);
        
        // Setting the active buffer to ping
        _isPingActive = 1;
    }

    /// <summary>
    /// Initializes the IDs of the ComputeShader's variables and initializes them.
    /// </summary>
    private void InitializeComputeShader()
    {
        // CONSTANTS
        
        // Simulation related constants
        simulationComputeShader.SetFloat("particles_amount", particlesAmount);
        simulationComputeShader.SetFloat("particle_mass", particleMass);
        
        // Caching the reference ID of the Shaders variables
        _halfBoundsID = Shader.PropertyToID("half_bounds");
        _collisionDampingID = Shader.PropertyToID("collision_damping");
        
        _gravityID = Shader.PropertyToID("gravity");

        _kernelRadiusID = Shader.PropertyToID("kernel_radius");
        _sqrKernelRadiusID = Shader.PropertyToID("sqr_kernel_radius");
        
        _poly6NormalizationID = Shader.PropertyToID("poly6_normalization");
        _spikyGradientFirstTermID = Shader.PropertyToID("spiky_gradient_first_term");
        _viscosityLaplacianFirstTermID = Shader.PropertyToID("viscosity_laplacian_first_term");

        _restDensityID = Shader.PropertyToID("rest_density");
        _stiffnessID = Shader.PropertyToID("stiffness");
        _viscosityID = Shader.PropertyToID("viscosity");

        _isPingActiveID = Shader.PropertyToID("is_ping_active");
        
        // Caching the reference ID of the ComputeShader's kernels
        _predictedPositionKernelID = simulationComputeShader.FindKernel("calculatePredictedPosition");
        _densityKernelID = simulationComputeShader.FindKernel("calculateDensity");
        _deltaVelocityKernelID = simulationComputeShader.FindKernel("calculateDeltaVelocity");
        
        // Assigning the buffer to the kernel
        simulationComputeShader.SetBuffer(_predictedPositionKernelID, "particles_buffer_ping", _particlesBufferPing);
        simulationComputeShader.SetBuffer(_predictedPositionKernelID, "particles_buffer_pong", _particlesBufferPong);
        simulationComputeShader.SetBuffer(_densityKernelID, "particles_buffer_ping", _particlesBufferPing);
        simulationComputeShader.SetBuffer(_densityKernelID, "particles_buffer_pong", _particlesBufferPong);
        simulationComputeShader.SetBuffer(_deltaVelocityKernelID, "particles_buffer_ping", _particlesBufferPing);
        simulationComputeShader.SetBuffer(_deltaVelocityKernelID, "particles_buffer_pong", _particlesBufferPong);

        UpdateComputeShaderVariables();
    }

    /// <summary>
    /// Initializes the ComputeShader's variables.
    /// </summary>
    private void UpdateComputeShaderVariables()
    {
        // Setting ComputeShader's vectors
        simulationComputeShader.SetVector(_halfBoundsID, _halfBounds);
        
        // Setting ComputeShader's floats
        simulationComputeShader.SetFloat(_collisionDampingID, collisionDamping);
        
        simulationComputeShader.SetFloat(_gravityID, gravity);
        
        simulationComputeShader.SetFloat(_kernelRadiusID, kernelRadius);
        simulationComputeShader.SetFloat(_sqrKernelRadiusID, _sqrKernelRadius);
        
        simulationComputeShader.SetFloat(_poly6NormalizationID, _poly6Normalization);
        simulationComputeShader.SetFloat(_spikyGradientFirstTermID, _spikyGradientNormalization);
        simulationComputeShader.SetFloat(_viscosityLaplacianFirstTermID, _viscosityLaplacianNormalization);
        
        simulationComputeShader.SetFloat(_restDensityID, restDensity);
        simulationComputeShader.SetFloat(_stiffnessID, stiffness);
        simulationComputeShader.SetFloat(_viscosityID, viscosity);
        
        simulationComputeShader.SetInt(_isPingActiveID, _isPingActive);
    }
    
    #endregion

    #region Unity Callback Functions

    /// <summary>
    /// Render the bounds of the simulation.
    /// </summary>
    private void OnDrawGizmos()
    {
        // Drawing the bounding gizmo
        Gizmos.color = Color.white;
        Gizmos.DrawWireCube(new Vector3(0, 0, 0), bounds);
    }

    /// <summary>
    /// Clean-up of the ComputeBuffers used for rendering
    /// </summary>
    private void OnDestroy()
    {
        // Releasing Compute Buffers to prevent memory leaks
        _particlesBufferPing?.Release();
        _particlesBufferPong?.Release();
    }

    /// <summary>
    /// Updates the underlying simulation.
    /// </summary>
    private void FixedUpdate()
    {
        // Updating the variables related to the kernels
        UpdateKernelVariables();
        
        // Updating the variables related to the bounds
        UpdateBoundsVariables();
        
        // Setting all the constants, in case they were changed via GUI
        UpdateComputeShaderVariables();
        
        // Computing a smaller delta time per substep
        var subDeltaTime = Time.fixedDeltaTime / simulationSubSteps;
        
        _gravityAcceleration = gravity * Vector2.down;
        
        // Running multiple simulation steps 
        for (var i = 0; i < simulationSubSteps; i++)
        {
            switch (integrationMethod)
            {
                default:
                case IntegrationMethod.Euler:
                    EulerSimulationStep(subDeltaTime);
                    break;
                case IntegrationMethod.LeapFrog:
                    LeapFrogSimulationStep(subDeltaTime);
                    break;
            }
        }

        // Returning 
        return;
        
        // TODO: TESTING GPU SIDE RENDERING
        
        // Computing the gravity vector
        // Running multiple simulation steps 
        for (var i = 0; i < simulationSubSteps; i++)
        {
            // Setting the delta time in the GPU
            simulationComputeShader.SetFloat("delta_time", subDeltaTime);
            
            // Dispatching the ComputeShader
            simulationComputeShader.Dispatch(_predictedPositionKernelID, particlesAmount / 256 + 1, 1, 1);
            simulationComputeShader.Dispatch(_densityKernelID, particlesAmount / 256 + 1, 1, 1);
            simulationComputeShader.Dispatch(_deltaVelocityKernelID, particlesAmount / 256 + 1, 1, 1);
            
            // Flipping the active buffer
            _isPingActive = _isPingActive == 1 ? 0 : 1;
            
            simulationComputeShader.SetInt(_isPingActiveID, _isPingActive);
        }
    }

    /// <summary>
    /// Renders the simulation.
    /// </summary>
    private void Update()
    {
        // If the simulation is not to be rendered, return
        if (!renderSimulation)
            return;
        
        // Updating the uniform Radius on GPU
        particleMaterial.SetFloat("particle_radius", particleRadius);
        particleMaterial.SetInt(_isPingActiveID, _isPingActive);

        // Updating position buffer on GPU
        // TODO: GPU SIDE RENDERING
        _particlesBufferPing.SetData(_particles);
        
        // Setting the rendering parameters
        var renderingParameters = new RenderParams(particleMaterial);
        
        // Render the particles
        Graphics.RenderPrimitives(renderingParameters, MeshTopology.Points, particlesAmount);
    }
    
    #endregion

    #region Simulation
    
    /// <summary>
    /// Execute a simulation step with the given delta time using Euler based integration.
    /// </summary>
    /// <param name="deltaTime">The delta time from the previous simulation step.</param>
    private void EulerSimulationStep(float deltaTime)
    {
        // Applying gravity and predict next positions
        Parallel.For(0, particlesAmount, currentIndex => 
        {
            // Applying gravity to the particle velocity
            _particles[currentIndex].velocity +=  _gravityAcceleration * deltaTime;
            
            // Predicting the next particle position
            _particles[currentIndex].predictedPosition = _particles[currentIndex].position 
                                                         + _particles[currentIndex].velocity * deltaTime;
        });
        
        // Updating the spatial lookup array
        UpdateSpatialLookup();
        
        // Calculating the density at particles position
        Parallel.For(0, particlesAmount, currentIndex =>
        {
            // Computing the density at particle position
            CalculateDensity(ref _particles[currentIndex]);
        });

        // Adding the current acceleration to the particle's velocity
        Parallel.For(0, particlesAmount, currentIndex =>
        {
            // Adding the acceleration to the particle velocity
            _particles[currentIndex].velocity += CalculateAcceleration(currentIndex) * deltaTime;
        });
        
        // Updating the particles position and resolving collision
        Parallel.For(0, particlesAmount, currentIndex =>
        {
            // Updating position
            _particles[currentIndex].position += _particles[currentIndex].velocity * deltaTime;

            // Resolving collisions
            ResolveCollisions(ref _particles[currentIndex].position, ref _particles[currentIndex].velocity);
        });
    }

    /// <summary>
    /// Execute a simulation step with the given delta time using Leapfrog based integration.
    /// </summary>
    /// <param name="deltaTime">The delta time from the previous simulation step.</param>
    private void LeapFrogSimulationStep(float deltaTime)
    {
        // Compute the half step velocity
        Parallel.For(0, particlesAmount, currentIndex => 
        {
            // Applying gravity to the particle velocity
            _particles[currentIndex].velocity +=  0.5f * deltaTime * _particles[currentIndex].acceleration;
        });
        
        // Computing the new position using the half step velocity
        Parallel.For(0, particlesAmount, currentIndex =>
        {
            // Updating position
            _particles[currentIndex].position += deltaTime * _particles[currentIndex].velocity;
        });
        
        // Updating the spatial lookup array
        UpdateSpatialLookup();
        
        // Calculating the density at particles position
        Parallel.For(0, particlesAmount, currentIndex =>
        {
            // Computing the density at particle position
            CalculateDensity(ref _particles[currentIndex]);
        });
        
        // Computing the acceleration at new position
        Parallel.For(0, particlesAmount, currentIndex =>
        {
            // Computing and the storing the acceleration at the half step position
            _particles[currentIndex].acceleration = CalculateAcceleration(currentIndex);
        });
        
        // Applying the acceleration to the velocity
        Parallel.For(0, particlesAmount, currentIndex =>
        {
            _particles[currentIndex].velocity += 0.5f * deltaTime * _particles[currentIndex].acceleration;
            
            // Resolving collisions
            ResolveCollisions(ref _particles[currentIndex].position, ref _particles[currentIndex].velocity);
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
        Parallel.For(0, particlesAmount, currentIndex => {
            // Mapping the particle position to a cell in the grid
            var cellCoordinates = integrationMethod == IntegrationMethod.Euler ?
                PositionToCellCoordinates(_particles[currentIndex].predictedPosition)
                : PositionToCellCoordinates(_particles[currentIndex].position);
            
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
        for (var i = 0; i < particlesAmount; i++)
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
        return Mathf.Abs(hash) % particlesAmount;
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
    /// <param name="targetParticle">The index of the particle to evaluate.</param>
    /// <returns>The density at the particle position</returns>
    private void CalculateDensity(ref Particle targetParticle)
    {
        // Computing the coordinates of the grid cell containing the given position
        var cellCoordinates = integrationMethod == IntegrationMethod.Euler? 
            PositionToCellCoordinates(targetParticle.predictedPosition)
            : PositionToCellCoordinates(targetParticle.position);
        
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
            for (var i = startIndex; i < particlesAmount && _spatialLookup[i].cellKey == currentCellKey; i++)
            {
                // Extracting a reference to the current particle
                ref var currentParticle = ref _particles[_spatialLookup[i].particleIndex];
                
                // Computing the distance between the current particle and the particle of interest
                var sqrDistance = 
                    integrationMethod == IntegrationMethod.Euler? 
                        (currentParticle.predictedPosition - targetParticle.predictedPosition).sqrMagnitude :
                        (currentParticle.position - targetParticle.position).sqrMagnitude;
            
                // Case in which the particle is outside the kernel radius or affecting itself
                if(sqrDistance >= _sqrKernelRadius)
                    continue;
            
                // Computing the influence of the current particle
                var influence = Poly6(sqrDistance);
            
                // Adding the result to the density
                density += particleMass * influence;
            }
        }
        
        targetParticle.density = density;
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
        // Extracting a reference to the current particle
        ref var targetParticle = ref _particles[targetParticleIndex];

        // Computing the coordinates of the grid cell containing the given position
        var cellCoordinates = integrationMethod == IntegrationMethod.Euler? 
            PositionToCellCoordinates(targetParticle.predictedPosition)
            : PositionToCellCoordinates(targetParticle.position);

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
            for (var i = startIndex; i < particlesAmount && _spatialLookup[i].cellKey == currentCellKey; i++)
            {
                // Extracting the index of the current particles
                var currentParticleIndex = _spatialLookup[i].particleIndex;
                
                // Skipping calculation for the particle by itself
                if (currentParticleIndex == targetParticleIndex)
                    continue;
                
                // Extracting a reference to the current particle
                ref var currentParticle = ref _particles[currentParticleIndex];
            
                // Computing the vector from target particle to current particle
                var targetToCurrent = integrationMethod == IntegrationMethod.Euler? 
                    targetParticle.predictedPosition - currentParticle.predictedPosition
                    : targetParticle.position - currentParticle.position;
                
                // Case in which the particle is outside the kernel radius
                if(targetToCurrent.sqrMagnitude >= _sqrKernelRadius)
                    continue;
                
                // Computing the distance
                var distance = targetToCurrent.magnitude;
                
                // Computing the influence of the current particle
                var viscosityInfluence = ViscosityLaplacian(distance);
                
                // Computing the difference between velocities
                var velocitiesDifference = currentParticle.velocity - targetParticle.velocity;
                
                // Adding the result to the force
                force += particleMass * viscosityInfluence * viscosity * velocitiesDifference 
                         / currentParticle.density;
                
                // Computing the direction between the current particle and the point of interest
                var pressureDirection = distance == 0 ?
                    new Vector3(1, 0, 0)
                    : targetToCurrent / distance;
            
                // Evaluating the slope of the kernel at the given distance
                var kernelSlope = SpikyGradient(distance);
            
                // Computing the shared pressure force
                var sharedPressure = CalculateSharedPressure(currentParticle.density, targetParticle.density);
                
                // Computing the applied force using both pressure and viscosity
                force -= particleMass * sharedPressure * kernelSlope * pressureDirection / currentParticle.density;
            }
        }
        return force / targetParticle.density + _gravityAcceleration;
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