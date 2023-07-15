//
// Created by Anchit Mishra on 2023-03-16.
//
#include <iostream>
#include <cmath>
#include <vector>
#include "Utils.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"

#ifndef FLUIDSIM_FLUID_H
#define FLUIDSIM_FLUID_H
#define LENGTH 1000
#define HEIGHT 1000
#define XDIM 60
#define YDIM 100
#define SPACING 10.0
#define PARTICLE_SIZE 0.25
#define SOLID 2
#define FLUID 0
#define EMPTY 1
#define PUSH_PENALTY 0.1
#define BETTER_PROJECTION true
#define BETTER_INTEGRATION true
#define TIMESTEP 1.0/60.0
#define SUBSTEPS 4
#define PENALTY false
#define PIC 0.1

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::ConjugateGradient;

class Fluid {
public:
    // Fluid should contain a grid of specified dimensions
    float spacing = SPACING; // This is the length of one cell in the grid
    int gridLength = floor(LENGTH / spacing) + 1;
    int gridHeight = floor(HEIGHT / spacing) + 1;
    int numCells = gridLength * gridHeight;

    // The grid is MAC - it stores velocities on the cube faces, pressure in the center
    // Current velocities
    std::vector<float> xVelocities;
    std::vector<float> yVelocities;
    // Marker array for extrapolation
    std::vector<int> xMarker;
    std::vector<int> yMarker;
    // Previous velocities
    std::vector<float> prevXVelocities;
    std::vector<float> prevYVelocities;
    // Normalization coefficients (for interpolation)
    std::vector<float> xR;
    std::vector<float> yR;
    // Each grid cell should also be marked as being solid, fluid or empty
    std::vector<int> cellType;

    // For Conjugate Gradient implementation
    float density = 1000;
    // To make sure drift doesn't occur (i.e. volume loss), we track particle density in each cell
    // This density must remain within a required range at all times, displacing any excess particles
    std::vector<float> cellParticleDensity;
    // We also need to track the initial particle density;
    float initialParticleDensity = 0.0; // We want roughly 4 particles per cell (2 along each dimension)
    // Spatial hashing table (unbounded grid)
    std::vector<int> spatialHashTable;
    std::vector<int> spatialHashParticles;
    int spatialHashTableSize;

    // Next up, we store penalty grid-related information for the fluid
    float particleRadius = spacing * PARTICLE_SIZE;
    float spatialHashGridSpacing = 2 * particleRadius;
    //int particleMassLength = floor((0.4 * LENGTH) / (2.0 * particleRadius));
    int particleMassLength = XDIM;
    //int particleMassHeight = floor((1.0 * HEIGHT) / (2.0 * particleRadius));
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

    // For surface reconstruction, we need to store grid point distances
    std::vector<float> pointDistance;
    // we also need to store indices of the closest point
    std::vector<int> closestPointIndices;

    // Finally, we come to functions
    Fluid();
    void simulateFluid();
    int IX(int x, int y);
    // Advection
    void advect();
    void detectBoundaryCollisions();
    void extrapolateVelocities();
    // Transfer
    void transferVelocitiesToGrid();
    void transferVelocitiesFromGrid();
    // Projection
    void projectGS();
    void projectPCG();
    void computeCellDensities();
    // Spatial Hashing (for collision detection)
    int spatialHashFunction(int x, int y);
    void spatialHashing();
    void detectParticleCollisions();
    // Surface reconstruction
    float distance(std::tuple<int, int> x, std::tuple<float, float> p);
    void computeSignedDistance();
};


#endif //FLUIDSIM_FLUID_H
