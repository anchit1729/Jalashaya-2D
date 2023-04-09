//
// Created by Anchit Mishra on 2023-04-07.
//
#include "Fluid.h"

// Hash function used as per description in 10 minute physics tutorial on spatial hashing
int Fluid::spatialHashFunction(int x, int y/*, int z*/) {
    return abs((x * 92837111) ^ (y * 689287499));// ^ (z * 283923481));
}

void Fluid::spatialHashing() {
    std::fill(spatialHashTable.begin(), spatialHashTable.end(), 0);
    // Compute which cell each particle falls into
    for (int i = 0; i < numParticles; i++)  {
        int cellXCoordinate = floor(particleXPositions[i] / spatialHashGridSpacing);
        int cellYCoordinate = floor(particleYPositions[i] / spatialHashGridSpacing);
        int cellIndex = spatialHashFunction(cellXCoordinate, cellYCoordinate) % (spatialHashTableSize - 1);
        // increase particle count at the specified index
        spatialHashTable[cellIndex]++;
    }

    // Iterate over the table and compute partial sums
    int accumulator = 0;
    for (int i = 0; i < spatialHashTableSize - 1; i++)  {
        accumulator += spatialHashTable[i];
        spatialHashTable[i] = accumulator;
    }
    // final entry stores the last accumulator value
    spatialHashTable[spatialHashTableSize - 1] = accumulator;

    // Now, populate the particle vector corresponding to the spatial hash
    for (int i = 0; i < numParticles; i++)  {
        int cellXCoordinate = floor(particleXPositions[i] / spatialHashGridSpacing);
        int cellYCoordinate = floor(particleYPositions[i] / spatialHashGridSpacing);
        int cellIndex = spatialHashFunction(cellXCoordinate, cellYCoordinate) % (spatialHashTableSize - 1);
        spatialHashTable[cellIndex]--;
        // store the index of the particle
        spatialHashParticles[spatialHashTable[cellIndex]] = i;
    }

}

void Fluid::detectParticleCollisions() {
    // of course, iterate over all particles
    for (int i = 0; i < numParticles; i++)  {
        // we want to detect particles within a distance of 2 * particleRadius
        int lowX = particleXPositions[i] - (2 * particleRadius); lowX = floor(lowX / spatialHashGridSpacing);
        int lowY = particleYPositions[i] - (2 * particleRadius); lowY = floor(lowY / spatialHashGridSpacing);
        int highX = particleXPositions[i] + (2 * particleRadius); highX = floor(highX / spatialHashGridSpacing);
        int highY = particleYPositions[i] + (2 * particleRadius); highY = floor(highY / spatialHashGridSpacing);
        for (int j = lowX; j < highX + 1; j++)  {
            for (int k = lowY; k < highY + 1; k++)   {
                int tableIndex = spatialHashFunction(j, k) % (spatialHashTableSize - 1);
                // check if particles are present in this index of the table
                int numParticlesInCell = spatialHashTable[tableIndex + 1] - spatialHashTable[tableIndex];
                if (numParticlesInCell > 0) {
                    // check particles
                    for (int a = spatialHashTable[tableIndex]; a < spatialHashTable[tableIndex] + numParticlesInCell; a++)  {
                        int particle = spatialHashParticles[a];
                        if (particle != i)  {
                            // compute distance
                            float x = (particleXPositions[i] - particleXPositions[particle]);
                            float y = (particleYPositions[i] - particleYPositions[particle]);
                            float distance = sqrt(pow(x, 2) + pow(y, 2));
                            if (distance < 2.0 * particleRadius && distance != 0)    {
                                // apply penalties, push particles (spring-like penalties)
                                particleXPositions[i] += PUSH_PENALTY * x;
                                particleYPositions[i] += PUSH_PENALTY * y;
                                particleXPositions[particle] -= PUSH_PENALTY * x;
                                particleYPositions[particle] -= PUSH_PENALTY * y;
                            }
                        }
                    }
                }
            }
        }
    }
}