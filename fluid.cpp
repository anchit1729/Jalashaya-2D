//
// Created by Anchit Mishra on 2023-03-16.
//

#include "fluid.h"

fluid::fluid() {
    // Part 1: Initialize grid stuff
    // Initialize the grid velocities
    xVelocities.resize(numCells, 0.0);
    yVelocities.resize(numCells, 0.0);
    // Initialize the marker values
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
            particleXPositions[index] = 25 * spacing + 2 * particleRadius * i; //+ (j % 2 == 0 ? particleRadius : 0);
            particleYPositions[index++] = 3 * spacing + 2 * particleRadius * j;
        }
    }

    containerWallXVelocity = 0;
    containerWallYVelocity = 0;

}

int fluid::IX(int x, int y) {
    return x * gridHeight + y;
}

float fluid::clamp(float val, float min, float max) {
    if (val < min) return min;
    if (val > max) return max;
    return val;
}

void fluid::simulateFluid() {
    // Implement simulation cycle here - transfer to grid, solve incompressibility, transfer from grid.
    detectBoundaryCollisions();
    transferVelocitiesToGrid();
    extrapolateVelocities();
    computeCellDensities();
    projectGS();
    extrapolateVelocities();
    transferVelocitiesFromGrid();
    advect(TIMESTEP/SUBSTEPS * 0.2);

}
