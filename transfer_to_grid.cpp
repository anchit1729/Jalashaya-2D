//
// Created by Anchit Mishra on 2023-04-07.
//
#include "fluid.h"

void fluid::transferVelocitiesToGrid() {
    // Transfer velocities of particles to their nearby grid side faces
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
            int cellIndex = IX(i, j);
            if (i == gridLength - 1 || j == gridHeight - 1 || i == 0 || j == 0) cellType[cellIndex] = SOLID;
            else cellType[cellIndex] = EMPTY;
        }
    }
    // Iterate over particles and mark fluid cells accordingly
    for (int i = 0; i < numParticles; i++)  {
        // Compute the cell coordinates
        int cellXCoordinate = clamp(floor(particleXPositions[i] / spacing), 0, gridLength - 1);
        int cellYCoordinate = clamp(floor(particleYPositions[i] / spacing), 0, gridHeight - 1);
        // Mark the cell as FLUID if not done already
        int cellIndex = cellXCoordinate * gridHeight + cellYCoordinate;

        if (cellXCoordinate < 0 || cellXCoordinate > gridLength - 1) std::cerr << "X out of bounds!\n";
        if (cellYCoordinate < 0 || cellYCoordinate > gridHeight - 1) std::cerr << "Y out of bounds!\n";
        if (cellType[cellIndex] == EMPTY) cellType[cellIndex] = FLUID;
    }
    // Now, commence transfer to grid
    // Part 1: Transfer X velocities
    float yShift = 0.5 * spacing; // Shift downwards to account for staggered grid
    for (int i = 0; i < numParticles; i++)  {
        // Again, compute the cell coordinates
        float x = clamp(particleXPositions[i], spacing, (gridLength - 1) * spacing);
        float y = clamp(particleYPositions[i], spacing, (gridHeight - 1) * spacing);
        int cellXCoordinate = floor(x / spacing);
        int cellYCoordinate = floor((y - yShift) / spacing);
        // the stuff below shouldn't happen
        if (cellXCoordinate < 0 || cellXCoordinate > gridLength - 1) std::cerr << "X out of bounds!\n";
        if (cellYCoordinate < 0 || cellYCoordinate > gridHeight - 1) std::cerr << "Y out of bounds!\n";
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
        if (corner1 < 0 || corner1 >= numCells) std::cerr << "1 Out of bounds!!\n";
        if (corner2 < 0 || corner2 >= numCells) std::cerr << "2 Out of bounds!!\n";
        if (corner3 < 0 || corner3 >= numCells) std::cerr << "3 Out of bounds!!\n";
        if (corner4 < 0 || corner4 >= numCells) std::cerr << "4 Out of bounds!!\n";
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
        // Finally, transfer particle velocity to the cell marks
        float particleVelocity = particleYVelocities[i];
        yVelocities[corner1] += particleVelocity * w1; yR[corner1] += w1; yMarker[corner1] = 0;
        yVelocities[corner2] += particleVelocity * w2; yR[corner2] += w2; yMarker[corner2] = 0;
        yVelocities[corner3] += particleVelocity * w3; yR[corner3] += w3; yMarker[corner3] = 0;
        yVelocities[corner4] += particleVelocity * w4; yR[corner4] += w4; yMarker[corner4] = 0;
    }
    // Wrap up the Y velocity transfer by restoring solid cell velocities as well as dividing by delta values (these are the summed weights from inverse bilinear interpolation)
    for (int i = 0; i < yVelocities.size(); i++)    {
        if (yR[i] > 0.0) {
            yVelocities[i] /= yR[i];
            yVelocities[i] += 9.81f * TIMESTEP/SUBSTEPS;
        }
    }
    //std::cout << "Velocity transfer complete.\n";

    // restore grid velocities for solid cells
    for (int i = 0; i < gridLength; i++) {
        for (int j = 0; j < gridHeight; j++) {
            if (cellType[IX(i, j)] == SOLID || cellType[IX(i - 1, j)] == SOLID)   {
                xVelocities[IX(i, j)] = prevXVelocities[IX(i, j)];
            }
            if (cellType[IX(i, j)] == SOLID || cellType[IX(i, j - 1)] == SOLID)   {
                yVelocities[IX(i, j)] = prevYVelocities[IX(i, j)];
            }
        }
    }
}