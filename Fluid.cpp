//
// Created by Anchit Mishra on 2023-03-16.
//

#include "Fluid.h"

Fluid::Fluid() {
    // Part 1: Initialize grid stuff
    // Initialize the grid velocities
    xVelocities.resize(numCells, 0.0);
    yVelocities.resize(numCells, 0.0);
    // Initialize the marker array
    xMarker.resize(numCells, INT_MAX);
    yMarker.resize(numCells, INT_MAX);
    // Initialize the previous grid velocities
    prevXVelocities.resize(numCells, 0.0);
    prevYVelocities.resize(numCells, 0.0);
    // Initialize the delta grid velocities
    xR.resize(numCells, 0.0);
    yR.resize(numCells, 0.0);
    // Initialize the cell type matrix
    cellType.resize(numCells, 0);
    // Initialize the cell particle density matrix
    cellParticleDensity.resize(numCells, 0.0);
    initialParticleDensity = 0;
    // Initialize the unbounded spatial hash table
    spatialHashTableSize = numParticles + 1;
    spatialHashTable.resize(2 * (spatialHashTableSize - 1) + 1, 0);
    spatialHashParticles.resize(numParticles);

    // For PCG Projection
    density = 1000.0;

    // Part 2: Initialize particle stuff
    // Initialize the particle positions
    particleXPositions.resize(numParticles, 0.0);
    particleYPositions.resize(numParticles, 0.0);
    // Initialize the particle velocities
    particleXVelocities.resize(numParticles, 0.0);
    particleYVelocities.resize(numParticles, 0.0);

    // Now that basic initialization is done, position the particles in a starting configuration for dam-break
    int index = 0;
    for (int i = 0; i < particleMassLength; i++)    {
        for (int j = 0; j < particleMassHeight; j++)    {
            particleXPositions[index] = spacing + 2 * particleRadius * i; //+ (j % 2 == 0 ? particleRadius : 0);
            particleYPositions[index++] = spacing + 2 * particleRadius * j;
        }
    }
    for (int i = 0; i < numParticles; i++) {
        float particleXCoordinate = 3 * spacing + 2 * particleRadius * (i % (gridLength)); // 2 particles per cell in X direction
        float particleYCoordinate = 3 * spacing + 2 * particleRadius * (i / (gridLength)); // 2 particles per cell in Y direction
        //particleXPositions[i] = particleXCoordinate;
        //particleYPositions[i] = particleYCoordinate;
    }
    containerWallXVelocity = 0;
    containerWallYVelocity = 0;

    pointDistance.resize(numCells, INT_MAX);
    closestPointIndices.resize(numCells, -1);
}

int Fluid::IX(int x, int y) {
    return x * gridHeight + y;
}

void Fluid::simulateFluid() {
    // Implement simulation cycle here - transfer to grid, solve incompressibility, transfer from grid.
    //std::cout << "Simulating fluid...\n";
    advect();
    detectBoundaryCollisions();
    if (PENALTY) {
        spatialHashing();
        detectParticleCollisions();
    }
    transferVelocitiesToGrid();
    extrapolateVelocities();
    if (BETTER_PROJECTION)  {
        projectPCG();
    } else {
        computeCellDensities();
        projectGS();
    }
    extrapolateVelocities();
    transferVelocitiesFromGrid();
    //std::cout << "Simulation step complete.\n";
}
