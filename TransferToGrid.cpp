//
// Created by Anchit Mishra on 2023-04-07.
//
#include "Fluid.h"

void Fluid::transferVelocitiesToGrid() {
    // Transfer velocities of particles to their nearby grid side faces
    //std::cout << "Transferring velocities to grid...\n";
    // First, as part of the FLIP method, we store a copy of all grid velocities
    prevXVelocities = xVelocities;
    prevYVelocities = yVelocities;
    // Now, reset the velocity values
    std::fill(xVelocities.begin(), xVelocities.end(), 0.0);
    std::fill(yVelocities.begin(), yVelocities.end(), 0.0);
    // Also reset the marker values
    std::fill(xMarker.begin(), xMarker.end(), INT_MAX);
    std::fill(yMarker.begin(), yMarker.end(), INT_MAX);
    // Also, make sure delta values are set to 0.0
    std::fill(xR.begin(), xR.end(), 0.0);
    std::fill(yR.begin(), yR.end(), 0.0);

    // Next, mark all boundary cells as solid cells to begin with
    for (int i = 0; i < gridLength; i++)  {
        for (int j = 0; j < gridHeight; j++)    {
            int cellIndex = i * gridHeight + j;
            if (i == 0 || i == gridLength - 1 || j == gridHeight - 1 || j == 0) cellType[cellIndex] = SOLID;
            else cellType[cellIndex] = EMPTY;
        }
    }
    // Iterate over particles and mark fluid cells accordingly
    for (int i = 0; i < numParticles; i++)  {
        // Compute the cell coordinates
        int cellXCoordinate = Utils::clamp(floor(particleXPositions[i] / spacing), 0, gridLength - 1);
        int cellYCoordinate = Utils::clamp(floor(particleYPositions[i] / spacing), 0, gridHeight - 1);
        // Mark the cell as FLUID if not done already
        int cellIndex = cellXCoordinate * gridHeight + cellYCoordinate;
        if (cellType[cellIndex] == EMPTY) cellType[cellIndex] = FLUID;
    }
    // Now, commence transfer to grid
    // Part 1: Transfer X velocities
    float yShift = 0.5 * spacing; // Shift downwards to account for staggered grid
    for (int i = 0; i < numParticles; i++)  {
        // Again, compute the cell coordinates
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
        // Finally, transfer particle velocity to the cell marks
        float particleVelocity = particleXVelocities[i];
        xVelocities[corner1] += particleVelocity * w1; xR[corner1] += w1; xMarker[corner1] = 0;
        xVelocities[corner2] += particleVelocity * w2; xR[corner2] += w2; xMarker[corner2] = 0;
        xVelocities[corner3] += particleVelocity * w3; xR[corner3] += w3; xMarker[corner3] = 0;
        xVelocities[corner4] += particleVelocity * w4; xR[corner4] += w4; xMarker[corner4] = 0;
    }
    // Wrap up the X velocity transfer by restoring solid cell velocities as well as dividing by delta values (these are the summed weights from inverse bilinear interpolation)
    for (int i = 0; i < xVelocities.size(); i++)    {
        if (xR[i] > 0.0) xVelocities[i] /= xR[i];
    }

    // Part 2: Transfer Y velocities
    float xShift = 0.5 * spacing; // Shift sideways to account for staggered grid
    for (int i = 0; i < numParticles; i++)  {
        // Again, compute the cell coordinates
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
        // Finally, transfer particle velocity to the cell marks
        float particleVelocity = particleYVelocities[i];
        yVelocities[corner1] += particleVelocity * w1; yR[corner1] += w1;
        yVelocities[corner2] += particleVelocity * w2; yR[corner2] += w2;
        yVelocities[corner3] += particleVelocity * w3; yR[corner3] += w3;
        yVelocities[corner4] += particleVelocity * w4; yR[corner4] += w4;
    }
    // Wrap up the Y velocity transfer by restoring solid cell velocities as well as dividing by delta values (these are the summed weights from inverse bilinear interpolation)
    for (int i = 0; i < yVelocities.size(); i++)    {
        if (yR[i] > 0.0) yVelocities[i] /= yR[i];
    }
    //std::cout << "Velocity transfer complete.\n";
}