//
// Created by Anchit Mishra on 2023-04-07.
//
#include "fluid.h"

void fluid::transferVelocitiesFromGrid() {
    // Transfer velocities from the grid faces to neighboring particles
    // First, transfer X velocities from the grid to particles
    float yShift = 0.5 * spacing;
    for (int i = 0; i < numParticles; i++)  {
        // Compute the cell coordinates
        float x = clamp(particleXPositions[i], spacing, (gridLength - 1) * spacing);
        float y = clamp(particleYPositions[i], spacing, (gridHeight - 1) * spacing);
        int cellXCoordinate = floor(x / spacing);
        int cellYCoordinate = floor((y - yShift) / spacing);
        // Extract the remainders (to compute corner weights w1, w2, w3, w4)
        float deltaX = x - cellXCoordinate * spacing;
        float deltaY = (y - yShift) - cellYCoordinate * spacing;
        // Now construct the corner weights
        float w1 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing));
        float w2 = (deltaX / spacing) * (1.0 - (deltaY / spacing));
        float w3 = (deltaX / spacing) * (deltaY / spacing);
        float w4 = (1.0 - (deltaX / spacing)) * (deltaY / spacing);
        // From here, compute the cell corners responsible for each weight
        int corner1 = cellXCoordinate * gridHeight + cellYCoordinate;
        int corner2 = fmin(cellXCoordinate + 1, gridLength - 2) * gridHeight + cellYCoordinate;
        int corner3 = fmin(cellXCoordinate + 1, gridLength - 2) * gridHeight + fmin(cellYCoordinate + 1, gridHeight - 2);
        int corner4 = cellXCoordinate * gridHeight + fmin(cellYCoordinate + 1, gridHeight - 2);
        // Proceed with transfer from the grid to particle
        // Only velocity values taken from SOLID or FLUID cells are valid
        float v1 = cellType[corner1] != EMPTY || cellType[corner1 - gridHeight] != EMPTY ? 1.0 : 0.0;
        float v2 = cellType[corner2] != EMPTY || cellType[corner2 - gridHeight] != EMPTY ? 1.0 : 0.0;
        float v3 = cellType[corner3] != EMPTY || cellType[corner3 - gridHeight] != EMPTY ? 1.0 : 0.0;
        float v4 = cellType[corner4] != EMPTY || cellType[corner4 - gridHeight] != EMPTY ? 1.0 : 0.0;
        float d = v1 * w1 + v2 * w2 + v3 * w3 + v4 * w4;
        float velocity = particleXVelocities[i];
        if (d > 0.0)    {
            float picVelocity = (v1 * w1 * xVelocities[corner1] + v2 * w2 * xVelocities[corner2] + v3 * w3 * xVelocities[corner3] + v4 * w4 * xVelocities[corner4]) / d;
            float weightedVelocityChanges = (v1 * w1 * (xVelocities[corner1] - prevXVelocities[corner1]) + v2 * w2 * (xVelocities[corner2] - prevXVelocities[corner2]) + v3 * w3 * (xVelocities[corner3] - prevXVelocities[corner3]) + v4 * w4 * (xVelocities[corner4] - prevXVelocities[corner4])) / d;
            float flipVelocity = velocity + weightedVelocityChanges;
            particleXVelocities[i] = PIC * picVelocity + (1 - PIC) * flipVelocity;
        }
    }

    // Next, transfer Y velocities from the grid to particles
    float xShift = 0.5 * spacing;
    for (int i = 0; i < numParticles; i++)  {
        // Compute the cell coordinates
        float x = clamp(particleXPositions[i], spacing, (gridLength - 1) * spacing);
        float y = clamp(particleYPositions[i], spacing, (gridHeight - 1) * spacing);
        int cellXCoordinate = floor((x - xShift) / spacing);
        int cellYCoordinate = floor(y / spacing);
        // Extract the remainders (to compute corner weights w1, w2, w3, w4)
        float deltaX = (x - xShift) - cellXCoordinate * spacing;
        float deltaY = y - cellYCoordinate * spacing;
        // Now construct the corner weights
        float w1 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing));
        float w2 = (deltaX / spacing) * (1.0 - (deltaY / spacing));
        float w3 = (deltaX / spacing) * (deltaY / spacing);
        float w4 = (1.0 - (deltaX / spacing)) * (deltaY / spacing);
        // From here, compute the cell corners responsible for each weight
        int corner1 = cellXCoordinate * gridHeight + cellYCoordinate;
        int corner2 = fmin(cellXCoordinate + 1, gridLength - 2) * gridHeight + cellYCoordinate;
        int corner3 = fmin(cellXCoordinate + 1, gridLength - 2) * gridHeight + fmin(cellYCoordinate + 1, gridHeight - 2);
        int corner4 = cellXCoordinate * gridHeight + fmin(cellYCoordinate + 1, gridHeight - 2);
        // Proceed with transfer from the grid to particle
        // Again, only SOLID or FLUID cell velocity values are valid for us
        float v1 = cellType[corner1] != EMPTY || cellType[corner1 - 1] != EMPTY ? 1.0 : 0.0;
        float v2 = cellType[corner2] != EMPTY || cellType[corner2 - 1] != EMPTY ? 1.0 : 0.0;
        float v3 = cellType[corner3] != EMPTY || cellType[corner3 - 1] != EMPTY ? 1.0 : 0.0;
        float v4 = cellType[corner4] != EMPTY || cellType[corner4 - 1] != EMPTY ? 1.0 : 0.0;
        float d = v1 * w1 + v2 * w2 + v3 * w3 + v4 * w4;
        float velocity = particleYVelocities[i];
        if (d > 0.0)    {
            float picVelocity = (v1 * w1 * yVelocities[corner1] + v2 * w2 * yVelocities[corner2] + v3 * w3 * yVelocities[corner3] + v4 * w4 * yVelocities[corner4]) / d;
            float weightedVelocityChanges = (v1 * w1 * (yVelocities[corner1] - prevYVelocities[corner1]) + v2 * w2 * (yVelocities[corner2] - prevYVelocities[corner2]) + v3 * w3 * (yVelocities[corner3] - prevYVelocities[corner3]) + v4 * w4 * (yVelocities[corner4] - prevYVelocities[corner4])) / d;
            float flipVelocity = velocity + weightedVelocityChanges;
            particleYVelocities[i] = PIC * picVelocity + (1 - PIC) * flipVelocity;
            if (x < 2 * spacing && particleYVelocities[i] < 0.1f) particleYVelocities[i] += 0.1f;
        }
    }
}