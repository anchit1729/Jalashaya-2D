//
// Created by Anchit Mishra on 2023-04-07.
//
#include "Fluid.h"

void Fluid::transferVelocitiesFromGrid() {
    // Transfer velocities from the grid faces to neighboring particles
    // First, transfer X velocities from the grid to particles
    float yShift = 0.5 * spacing;
    for (int i = 0; i < numParticles; i++)  {
        // Compute the cell coordinates
        float x = Utils::clamp(particleXPositions[i], spacing, (gridLength - 1) * spacing);
        float y = Utils::clamp(particleYPositions[i], spacing, (gridHeight - 1) * spacing);
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
        int offset = gridHeight;
        // Determine whether contributing corners on grid are valid or not
        float valid1 = cellType[corner1] != EMPTY || cellType[corner1 - offset] != EMPTY ? 1.0 : 0.0;
        float valid2 = cellType[corner2] != EMPTY || cellType[corner2 - offset] != EMPTY ? 1.0 : 0.0;
        float valid3 = cellType[corner3] != EMPTY || cellType[corner3 - offset] != EMPTY ? 1.0 : 0.0;
        float valid4 = cellType[corner4] != EMPTY || cellType[corner4 - offset] != EMPTY ? 1.0 : 0.0;
        float d = w1 * valid1 + w2 * valid2 + w3 * valid3 + w4 * valid4;
        float velocity = particleXVelocities[i];
        if (d > 0.0)    {
            float picVelocity = (valid1 * w1 * xVelocities[corner1] + valid2 * w2 * xVelocities[corner2] + valid3 * w3 * xVelocities[corner3] + valid4 * w4 * xVelocities[corner4]) / d;
            float weightedVelocityChanges = (valid1 * w1 * (xVelocities[corner1] - prevXVelocities[corner1]) + valid2 * w2 * (xVelocities[corner2] - prevXVelocities[corner2]) + valid3 * w3 * (xVelocities[corner3] - prevXVelocities[corner3]) + valid4 * w4 * (xVelocities[corner4] - prevXVelocities[corner4])) / d;
            float flipVelocity = velocity + weightedVelocityChanges;
            particleXVelocities[i] = PIC * picVelocity + (1 - PIC) * flipVelocity;
        }
    }

    // Next, transfer Y velocities from the grid to particles
    float xShift = 0.5 * spacing;
    for (int i = 0; i < numParticles; i++)  {
        // Compute the cell coordinates
        float x = Utils::clamp(particleXPositions[i], spacing, (gridLength - 1) * spacing);
        float y = Utils::clamp(particleYPositions[i], spacing, (gridHeight - 1) * spacing);
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
        int offset = 1;
        // Determine whether contributing corners on grid are valid or not
        float valid1 = cellType[corner1] != EMPTY || cellType[corner1 - offset] != EMPTY ? 1.0 : 0.0;
        float valid2 = cellType[corner2] != EMPTY || cellType[corner2 - offset] != EMPTY ? 1.0 : 0.0;
        float valid3 = cellType[corner3] != EMPTY || cellType[corner3 - offset] != EMPTY ? 1.0 : 0.0;
        float valid4 = cellType[corner4] != EMPTY || cellType[corner4 - offset] != EMPTY ? 1.0 : 0.0;
        float d = w1 * valid1 + w2 * valid2 + w3 * valid3 + w4 * valid4;
        float velocity = particleYVelocities[i];
        if (d > 0.0)    {
            float picVelocity = (valid1 * w1 * yVelocities[corner1] + valid2 * w2 * yVelocities[corner2] + valid3 * w3 * yVelocities[corner3] + valid4 * w4 * yVelocities[corner4]) / d;
            float weightedVelocityChanges = (valid1 * w1 * (yVelocities[corner1] - prevYVelocities[corner1]) + valid2 * w2 * (yVelocities[corner2] - prevYVelocities[corner2]) + valid3 * w3 * (yVelocities[corner3] - prevYVelocities[corner3]) + valid4 * w4 * (yVelocities[corner4] - prevYVelocities[corner4])) / d;
            float flipVelocity = velocity + weightedVelocityChanges;
            particleYVelocities[i] = PIC * picVelocity + (1 - PIC) * flipVelocity;
        }
    }
}