# Efficient Fluid Simulation Using Smoothed Particle Hydrodynamics in Unity

Real-time, GPU-accelerated SPH fluids in Unity 6 with volumetric rendering (URP), designed and evaluated as my Bachelor’s thesis at [USI](https://www.usi.ch/en).

![Particles](https://github.com/GuglielmoMazzesiDaniele/SPH/blob/main/Images/Particles.png)

## TL;DR

- **Method**: Weakly-compressible SPH with density/pressure, viscosity, near-density stabilization, and symplectic Euler integration. All core steps run in Compute Shaders.
- **Acceleration**: GPU spatial hashing + prefix sums for O(n) neighbor queries (3×3×3 cells).
- **Rendering**: URP ray-marched volume with Beer’s Law absorption, Schlick Fresnel reflection, and normals from the density gradient.
- **Performance** (test rig below): ≈ 600k particles ~ 60 FPS; up to 1.2M particles near 30 FPS.

## Demo

![Demo](https://github.com/GuglielmoMazzesiDaniele/SPH/blob/main/GIFs/Basic%20Showcase.gif)

## Features
- **End-to-end GPU pipeline**: prediction → spatial grid build/sort → density/pressure/viscosity → integration.
- **Stable field evaluation** using predicted positions + near-pressure term to reduce clumping and induce surface-tension-like behavior.
- **Volumetric look**: single cube bounds with ray marching through a 3D density field; depth-aware absorption + Fresnel blending.
- **Interactive controls** for pausing, resetting, screenshots, slow-mo, impulses, and free-fly camera.

## How it works

### SPH Solver (physics)
- **Density/pressure**: symmetric SPH formulation with an equation of state; viscosity via velocity differences and kernel Laplacian; gravity as body force. Integrated using symplectic Euler.
- **Stabilization**:
  - *Predicted positions* for field sampling (reduces feedback explosions).
  - *Near-pressure* short-range repulsion to prevent clumping and add cohesion at the free surface.

### GPU pipeline (performance)
- **Uniform grid + spatial hashing** built in compute: hash, histogram, exclusive prefix sum, scatter/sort, then per-cell offsets for constant-time lookups.
- **Hashing strategy**: a block-partitioned linear combination of cell coords minimized collisions and outperformed modular/Morton variants in my benchmarks.
 
### Rendering (URP)
- **Ray marching** inside a bounding cube, sampling a 3D density texture updated each frame. Front-to-back compositing with early-exit on opacity.
- **Beer’s Law** (density-weighted): transmittance to model depth/opacity.
- **Fresnel (Schlick's Approximation)**: reflection vs. refraction weighting by view angle; reflections at glancing angles, transparency head-on.
- **Normals**: central-difference gradient of density for lighting/Fresnel.

![Honey](https://github.com/GuglielmoMazzesiDaniele/SPH/blob/main/GIFs/Honey.gif)

## Requirements
- Unity 6 with Universal Render Pipeline (URP).
- GPU with Compute Shader support (DX11+/Metal/Vulkan). Tested on GeForce RTX 3070 (8 GB).
- Windows 10 used for benchmarks; other OSes should work if compute is supported.

## Project Setup
- Clone the repo and open it in Unity 6 (URP project).
- Open the sample scene (e.g., FluidDemo.unity).
- Press Play. You can use the controls below to interact.
- Optional: tweak Simulation and Rendering parameters in the inspector (kernel radius, stiffness, substeps, volume resolution, absorption, Fresnel base reflectivity). Defaults are listed under Configuration.

## Controls
- Space – Play/Pause
- R – Reset
- X / Z – Single / burst screenshots
- F – Slow-motion toggle
- Arrow Keys – Apply directional force to the fluid
- W/A/S/D, Q/E – Fly camera (shift for turbo)
- Mouse Right-Drag look, Middle-Drag pan, Scroll zoom
*See the thesis for the full mapping*.



