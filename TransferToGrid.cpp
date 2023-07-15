//
// Created by Anchit Mishra on 2023-04-07.
//
#include "Fluid.h"

void Fluid::transferVelocitiesToGrid() {
    // Transfer velocities of particles to their nearby grid side faces
    // First, store a copy of previous velocity field values (for FLIP)
    prevXVelocities = xVelocities;
    prevYVelocities = yVelocities;
    std::fill(xVelocities.begin(), xVelocities.end(), 0);
    std::fill(yVelocities.begin(), yVelocities.end(), 0);
    std::fill(xR.begin(), xR.end(), 0);
    std::fill(yR.begin(), yR.end(), 0);

    // Next, mark all boundary cells as solid cells
    for (int i = 0; i < gridLength; i++)  {
        for (int j = 0; j < gridHeight; j++)    {
            int cellIX = IX(i, j);
            if (i == 0 || i == gridLength - 1 || j == 0 || j == gridHeight - 1) cellType[cellIX] = SOLID;
            else cellType[cellIX] = EMPTY;
        }
    }

    // Mark all EMPTY cells that contain fluid
    for (int i = 0; i < numParticles; i++)  {
        int cellXCoordinate = Utils::clamp(floor(particleXPositions[i] / spacing), 0, gridLength - 1);
        int cellYCoordinate = Utils::clamp(floor(particleYPositions[i] / spacing), 0, gridHeight - 1);
        int cellIX = IX(cellXCoordinate, cellYCoordinate);
        if (cellType[cellIX] == EMPTY) cellType[cellIX] = FLUID;
    }

    // Now, begin with transferring all x velocity components
    float yShift = 0.5 * spacing; // Bring the grid up/down by half the spacing to align corners with cell centres
    for (int i = 0; i < numParticles; i++)  {
        float x = Utils::clamp(particleXPositions[i], spacing, (gridLength - 1) * spacing);
        float y = Utils::clamp(particleYPositions[i], spacing, (gridHeight - 1) * spacing);
        int cellXCoordinate = floor(x / spacing);
        int cellYCoordinate = floor((y - yShift) / spacing);
        // Extract the remainders (to compute corner weights w1, w2, w3, w4)
        float deltaX = x - cellXCoordinate * spacing;
        float deltaY = (y - yShift) - cellYCoordinate * spacing;
        // Compute corner weights
        float w1 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing));
        float w2 = (deltaX / spacing) * (1.0 - (deltaY / spacing));
        float w3 = (deltaX / spacing) * (deltaY / spacing);
        float w4 = (1.0 - (deltaX / spacing)) * (deltaY / spacing);
        // Next, compute the cell corners
        int c1 = IX(cellXCoordinate, cellYCoordinate);
        int c2 = IX(cellXCoordinate + 1, cellYCoordinate);
        int c3 = IX(cellXCoordinate + 1, cellYCoordinate + 1);
        int c4 = IX(cellXCoordinate, cellYCoordinate + 1);
        // Finally perform velocity transfer
        float v = particleXVelocities[i];
        xVelocities[c1] = v * w1; xR[c1] += w1;
        xVelocities[c2] = v * w2; xR[c2] += w2;
        xVelocities[c3] = v * w3; xR[c3] += w3;
        xVelocities[c4] = v * w4; xR[c4] += w4;
    }

    for (int i = 0; i < numCells; i++)  {
        if (xR[i] > 0) xVelocities[i] /= xR[i];
    }

    // Next up, transfer y velocity components
    float xShift = 0.5 * spacing;
    for (int i = 0; i < numParticles; i++)  {
        float x = Utils::clamp(particleXPositions[i], spacing, (gridLength - 1) * spacing);
        float y = Utils::clamp(particleYPositions[i], spacing, (gridHeight - 1) * spacing);
        int cellXCoordinate = floor((x - xShift) / spacing);
        int cellYCoordinate = floor(y / spacing);
        // Extract the remainders (to compute corner weights w1, w2, w3, w4)
        float deltaX = (x - xShift) - cellXCoordinate * spacing;
        float deltaY = y - cellYCoordinate * spacing;
        // Compute corner weights
        float w1 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing));
        float w2 = (deltaX / spacing) * (1.0 - (deltaY / spacing));
        float w3 = (deltaX / spacing) * (deltaY / spacing);
        float w4 = (1.0 - (deltaX / spacing)) * (deltaY / spacing);
        // Next, compute cell corners
        int c1 = IX(cellXCoordinate, cellYCoordinate);
        int c2 = IX(cellXCoordinate + 1, cellYCoordinate);
        int c3 = IX(cellXCoordinate + 1, cellYCoordinate + 1);
        int c4 = IX(cellXCoordinate, cellYCoordinate + 1);
        // Finally, perform velocity transfer
        float v = particleYVelocities[i];
        yVelocities[c1] = v * w1; yR[c1] += w1;
        yVelocities[c2] = v * w2; yR[c2] += w2;
        yVelocities[c3] = v * w3; yR[c3] += w3;
        yVelocities[c4] = v * w4; yR[c4] += w4;
    }

    for (int i = 0; i < numCells; i++)  {
        if (yR[i] > 0) yVelocities[i] /= yR[i];
    }

    // Finally, restore solid cell velocities to zero
    for (int i = 0; i < gridLength; i++)    {
        for (int j = 0; j < gridHeight; j++)    {
            if (cellType[IX(i, j)] == SOLID || cellType[IX(i - 1, j)] == SOLID) xVelocities[IX(i, j)] = prevXVelocities[IX(i, j)];
            if (cellType[IX(i, j)] == SOLID || cellType[IX(i, j - 1)] == SOLID) yVelocities[IX(i, j)] = prevYVelocities[IX(i, j)];
        }
    }

}