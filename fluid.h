//
// Created by Anchit Mishra on 2023-03-16.
//
#include <iostream>
#include <cmath>
#include <vector>


#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"


#ifndef FLUIDSIM_FLUID_H
#define FLUIDSIM_FLUID_H
// Defining grid dimensions
#define LENGTH 800
#define HEIGHT 800
// Defining the number of fluid particles along each dimension
#define XDIM 100
#define YDIM 100
// Defining the grid cell spacing
#define SPACING 10.0
// Defining the multiplication factor for particle size in terms of grid cell spacing (each particle is PARTICLE_SIZE * SPACING in radius)
#define PARTICLE_SIZE 0.22
#define FLUID 0
#define EMPTY 1
#define SOLID 2
// Defining whether to use PCG (better but slower but also parallel) or GS (standard)
#define BETTER_PROJECTION true
// Defining the time-step size
#define TIMESTEP 1.0/60.0
// Defining the number of sub-steps (to satisfy CFL condition)
#define SUBSTEPS 1
// Defining the PIC/FLIP blending ratio
#define PUSH_PENALTY 0.1
#define BETTER_INTEGRATION true
#define TIMESTEP 1.0/60.0
#define SUBSTEPS 0.5
#define PENALTY false
#define PIC 0.05

class fluid {
public:
    // fluid should contain a grid of specified dimensions
    float spacing = SPACING; // This is the length of one cell in the grid
    int gridLength = floor(LENGTH / spacing) + 1; // The number of grid cells along the X direction
    int gridHeight = floor(HEIGHT / spacing) + 1; // The number of grid cells along the Y direction
    int numCells = gridLength * gridHeight;

    // The grid is MAC - it stores velocities on the cube faces, pressure in the center
    // Pressure values are computed on the fly as part of projection
    // Current velocities
    std::vector<float> xVelocities; // Stores the X velocity of the grid cell face on the left
    std::vector<float> yVelocities; // Stores the Y velocity of the grid cell face on the right
    // Marker values
    std::vector<int> xMarker;
    std::vector<int> yMarker;
    // Previous velocities
    std::vector<float> prevXVelocities; // Similar to regular xVelocities vector
    std::vector<float> prevYVelocities; // Similar to regular yVelocities vector
    // Normalization coefficients (for interpolation)
    std::vector<float> xR;
    std::vector<float> yR;
    // Each grid cell should also be marked as being solid, fluid or empty
    std::vector<int> cellType;

    // For Conjugate Gradient implementation
    float density = 1000; // Use true-to-life values for water
    // To make sure drift doesn't occur (i.e. volume loss), we track particle density in each cell
    // This density must remain within a required range at all times, displacing any excess particles
    std::vector<float> cellParticleDensity;
    // We also need to track the initial particle density;
    float initialParticleDensity = 0.0; // We want roughly 4 particles per cell (2 along each dimension)

    // Next up, we store penalty grid-related information for the fluid
    float particleRadius = spacing * PARTICLE_SIZE;
    int particleMassLength = XDIM;
    int particleMassHeight = YDIM;
    int numParticles = particleMassHeight * particleMassLength;
    // Store particle positions
    std::vector<float> particleXPositions;
    std::vector<float> particleYPositions;
    // Store particle velocities
    std::vector<float> particleXVelocities;
    std::vector<float> particleYVelocities;
    float containerWallXVelocity;
    float containerWallYVelocity;

    // Finally, we come to functions
    fluid();
    void simulateFluid();
    int IX(int x, int y);
    // Advection
    void advect(float dt);
    void extrapolateVelocities();
    void detectBoundaryCollisions();
    // Transfer
    void transferVelocitiesToGrid();
    void transferVelocitiesFromGrid();
    // Projection
    void projectGS();
    void computeCellDensities();
    // Utility function for clamping
    static float clamp(float val, float min, float max);
};


#endif //FLUIDSIM_FLUID_H
