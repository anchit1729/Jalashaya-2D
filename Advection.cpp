//
// Created by Anchit Mishra on 2023-04-07.
//
#include "Fluid.h"

void Fluid::advect() {
    float dt = TIMESTEP/SUBSTEPS;
    for (int i = 0; i < numParticles; i++)  {
        // Only velocity in the Y direction is integrated - due to gravity
        particleYVelocities[i] += dt * 9.81f;
        particleXPositions[i] += dt * particleXVelocities[i];
        particleYPositions[i] += dt * particleYVelocities[i];
    }
}

void Fluid::detectBoundaryCollisions() {
    // First, set a minimum and maximum limit for particles to be located at
    float minX = 2 * particleRadius;
    float minY = 2 * particleRadius;
    float maxX = spacing * (gridLength - 1) - particleRadius;
    float maxY = spacing * (gridHeight - 1) - particleRadius;
    // Iterate over all particles, see whether they are colliding with any boundaries
    for (int i = 0; i < numParticles; i++)  {
        // Clamp the X and Y coordinates of the particle between boundary coordinates
        if (particleXPositions[i] < minX || particleXPositions[i] > maxX) particleXVelocities[i] = containerWallXVelocity;
        if (particleYPositions[i] < minY || particleYPositions[i] > maxY) particleYVelocities[i] = containerWallYVelocity;
        particleXPositions[i] = Utils::clamp(particleXPositions[i], minX, maxX);
        particleYPositions[i] = Utils::clamp(particleYPositions[i], minY, maxY);
    }
}
