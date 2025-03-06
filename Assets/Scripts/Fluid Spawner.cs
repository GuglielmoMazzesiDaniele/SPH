using System;
using System.Collections.Generic;
using UnityEngine;

public class FluidSpawner : MonoBehaviour
{
    [Header("Spawn Regions")]
    [SerializeField] public SpawnRegion[] spawnRegions;

    #region Unity Callback Functions

    private void OnDrawGizmos()
    {
        // If the simulation is running, return
        if (Application.isPlaying)
            return;
        
        // Iterating all spawn regions
        foreach (var spawnRegion in spawnRegions)
        {
            // Drawing the spawn region bounding box
            Gizmos.color = spawnRegion.boundariesColor;
            Gizmos.DrawWireCube(spawnRegion.position, spawnRegion.size);
        }
    }
    
    #endregion
    
    /// <summary>
    /// Returns an array containing all the particles spawned by this spawner.
    /// </summary>
    /// <returns>A Vector3 array containing all the starting particles positions.</returns>
    public Vector3[] GetSpawnPositions()
    {
        // Initializing the total amount of positions
        var totalPositions = new List<Vector3>();
        
        // Iterating all the spawn points and generating particles using their data
        foreach (var spawnRegion in spawnRegions)
        {
            // Adding the particle positions of the spawn region to the total amount of positions
            totalPositions.AddRange(spawnRegion.ParticlePositions());
        }

        // Returning the list as an array
        return totalPositions.ToArray();
    }

    /// <summary>
    /// Returns the total amount of particles spawned by the spawner.
    /// </summary>
    /// <returns>An integer representing the total amount of particles</returns>
    public int GetParticlesAmount()
    {
        // Initializing the particles amount
        var particlesAmount = 0;
        
        // Iterating all spawn regions
        foreach (var spawnRegion in spawnRegions)
        {
            // Getting the particles per axis
            var particlesPerAxis = spawnRegion.ParticlesPerAxis();
            
            // Increasing the particlesAmount
            particlesAmount += particlesPerAxis.x * particlesPerAxis.y * particlesPerAxis.z;
        }

        return particlesAmount;
    }
    
    // Struct representing a spawn region
    [Serializable]
    public class SpawnRegion
    {
        public uint particlesAmount;
        public Vector3 position;
        public Vector3 size;
        public Color boundariesColor;

        /// <summary>
        /// Returns the list of all the particles belonging to this spawn region.
        /// </summary>
        /// <returns>A List of Vector3 containing the positions.</returns>
        public List<Vector3> ParticlePositions()
        {
            // Initializing the list of particles
            var positions = new List<Vector3>();
            
            // Computing the particles per axis
            var particlesPerAxis = ParticlesPerAxis();
            
            // Filling the particles list
            for (var currentX = 0; currentX < particlesPerAxis.x; currentX++)
            {
                for (var currentY = 0; currentY < particlesPerAxis.y; currentY++)
                {
                    for (var currentZ = 0; currentZ < particlesPerAxis.z; currentZ++)
                    {
                        // Computing the X coordinate
                        var positionX = (currentX / Mathf.Max(particlesPerAxis.x - 1f, 1.0f) - 0.5f) * size.x 
                                        + position.x;
                        // Computing the Y coordinate
                        var positionY = (currentY / Mathf.Max(particlesPerAxis.y - 1f, 1.0f) - 0.5f) * size.y 
                                        + position.y;
                        // Computing the Z coordinate
                        var positionZ = (currentZ / Mathf.Max(particlesPerAxis.z - 1f, 1.0f) - 0.5f) * size.z 
                                        + position.z;
                        
                        // Pushing the position on the list
                        positions.Add(new Vector3(positionX, positionY, positionZ));
                    }
                }
            }

            return positions;
        }

        /// <summary>
        /// Generates a Vector3Int containing the amount of particles per axis as close as possible to
        /// the desired particles amount.
        /// </summary>
        /// <returns>A Vector3Int containing the amount of particles per axis.</returns>
        public Vector3Int ParticlesPerAxis()
        {
            if (particlesAmount == 0)
                return new Vector3Int(0, 0, 0);
            
            // Finding the biggest axis of the spawn region
            var maxSize = Mathf.Max(size.x, size.y, size.z);
            
            // Computing the ratios of particles per axis with respect to the maximum
            var ratioPerAxis = new Vector3(
                size.x / maxSize, 
                size.y / maxSize, 
                size.z / maxSize);
            
            // Auxiliary number used to compute the amount per axis
            var startingAmount = Mathf.Pow(particlesAmount, 1f / 3f);

            // Computing the ideal number of particles per axis
            var amountPerAxis = new Vector3Int(
                Mathf.Max(1, Mathf.RoundToInt(startingAmount * ratioPerAxis.x)),
                Mathf.Max(1, Mathf.RoundToInt(startingAmount * ratioPerAxis.y)),
                Mathf.Max(1, Mathf.RoundToInt(startingAmount * ratioPerAxis.z))
            );
            
            // Computing the current amount of total particles
            var computedTotal = amountPerAxis.x * amountPerAxis.y * amountPerAxis.z;
            
            // Reducing the amount to get it as close at possible to the desired amount
            while (computedTotal > particlesAmount)
            {
                // Removing particles 
                if (amountPerAxis.x > 1) amountPerAxis.x--;
                else if (amountPerAxis.y > 1) amountPerAxis.y--;
                else if (amountPerAxis.z > 1) amountPerAxis.z--;
                
                // Recomputing the total
                computedTotal = amountPerAxis.x * amountPerAxis.y * amountPerAxis.z;
            }
            
            // Increasing the amount to get it as close at possible to the desired amount
            while (computedTotal < particlesAmount)
            {
                // Adding particles 
                if (amountPerAxis.x <= amountPerAxis.y && amountPerAxis.x <= amountPerAxis.z) amountPerAxis.x++;
                else if (amountPerAxis.y <= amountPerAxis.x && amountPerAxis.y <= amountPerAxis.z) amountPerAxis.y++;
                else amountPerAxis.z++;
                
                // Recomputing the total
                computedTotal = amountPerAxis.x * amountPerAxis.y * amountPerAxis.z;
            }

            return amountPerAxis;
        }
    }
}