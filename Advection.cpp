//
// Created by Anchit Mishra on 2023-04-07.
//
#include "Fluid.h"

void Fluid::advect() {
    float dt = TIMESTEP/SUBSTEPS;
    // First, integrate the velocities and then integrate the positions
    //We try out an RK3 setup
    if (BETTER_INTEGRATION) {
        for (int i = 0; i < numParticles; i++)  {
            // Extract stage 1 values, k1_i
            float k1_x = particleXVelocities[i];
            float k1_y = particleYVelocities[i];
            // Now comes the hard part - derive stage 2 velocities
            float k2_x = 0;
            float k2_y = 0;
            float secondXPosition = particleXPositions[i] + 0.5 * dt * k1_x;
            float secondYPosition = particleYPositions[i] + 0.5 * dt * k1_y;
            float yShift = 0.5 * spacing;
            // Interpolate the velocity field to get velocities at this position
            int cellXCoordinate = floor(secondXPosition / spacing);
            int cellYCoordinate = floor((secondYPosition - yShift) / spacing);
            // Extract the remainders (to compute corner weights w1, w2, w3, w4)
            float deltaX = secondXPosition - cellXCoordinate * spacing;
            float deltaY = (secondYPosition - yShift) - cellYCoordinate * spacing;
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
            float valid1 = cellType[corner1] != EMPTY || cellType[corner1 - offset] != EMPTY ? 1.0 : 1.0;
            float valid2 = cellType[corner2] != EMPTY || cellType[corner2 - offset] != EMPTY ? 1.0 : 1.0;
            float valid3 = cellType[corner3] != EMPTY || cellType[corner3 - offset] != EMPTY ? 1.0 : 1.0;
            float valid4 = cellType[corner4] != EMPTY || cellType[corner4 - offset] != EMPTY ? 1.0 : 1.0;
            float d = w1 * valid1 + w2 * valid2 + w3 * valid3 + w4 * valid4;
            float velocity = particleXVelocities[i];
            if (d > 0.0)    {
                float picVelocity = (valid1 * w1 * xVelocities[corner1] + valid2 * w2 * xVelocities[corner2] + valid3 * w3 * xVelocities[corner3] + valid4 * w4 * xVelocities[corner4]) / d;
                float weightedVelocityChanges = (valid1 * w1 * (xVelocities[corner1] - prevXVelocities[corner1]) + valid2 * w2 * (xVelocities[corner2] - prevXVelocities[corner2]) + valid3 * w3 * (xVelocities[corner3] - prevXVelocities[corner3]) + valid4 * w4 * (xVelocities[corner4] - prevXVelocities[corner4])) / d;
                float flipVelocity = velocity + weightedVelocityChanges;
                k2_x = PIC * picVelocity + (1 - PIC) * flipVelocity;
            }

            float xShift = 0.5 * spacing;
            // Do the same for stage 2 y velocities
            cellXCoordinate = floor((secondXPosition - xShift) / spacing);
            cellYCoordinate = floor(secondYPosition / spacing);
            // Extract the remainders (to compute corner weights w1, w2, w3, w4)
            deltaX = (secondXPosition - xShift) - cellXCoordinate * spacing;
            deltaY = secondYPosition - cellYCoordinate * spacing;
            // Now construct the corner weights
            w1 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing));
            w2 = (deltaX / spacing) * (1.0 - (deltaY / spacing));
            w3 = (deltaX / spacing) * (deltaY / spacing);
            w4 = (1.0 - (deltaX / spacing)) * (deltaY / spacing);
            // From here, compute the cell corners responsible for each weight
            corner1 = cellXCoordinate * gridHeight + cellYCoordinate;
            corner2 = fmin(cellXCoordinate + 1, gridLength - 2) * gridHeight + cellYCoordinate;
            corner3 = fmin(cellXCoordinate + 1, gridLength - 2) * gridHeight + fmin(cellYCoordinate + 1, gridHeight - 2);
            corner4 = cellXCoordinate * gridHeight + fmin(cellYCoordinate + 1, gridHeight - 2);
            // Proceed with transfer from the grid to particle
            offset = 1;
            // Determine whether contributing corners on grid are valid or not
            valid1 = cellType[corner1] != EMPTY || cellType[corner1 - offset] != EMPTY ? 1.0 : 1.0;
            valid2 = cellType[corner2] != EMPTY || cellType[corner2 - offset] != EMPTY ? 1.0 : 1.0;
            valid3 = cellType[corner3] != EMPTY || cellType[corner3 - offset] != EMPTY ? 1.0 : 1.0;
            valid4 = cellType[corner4] != EMPTY || cellType[corner4 - offset] != EMPTY ? 1.0 : 1.0;
            d = w1 * valid1 + w2 * valid2 + w3 * valid3 + w4 * valid4;
            velocity = particleYVelocities[i];
            if (d > 0.0)    {
                float picVelocity = (valid1 * w1 * yVelocities[corner1] + valid2 * w2 * yVelocities[corner2] + valid3 * w3 * yVelocities[corner3] + valid4 * w4 * yVelocities[corner4]) / d;
                float weightedVelocityChanges = (valid1 * w1 * (yVelocities[corner1] - prevYVelocities[corner1]) + valid2 * w2 * (yVelocities[corner2] - prevYVelocities[corner2]) + valid3 * w3 * (yVelocities[corner3] - prevYVelocities[corner3]) + valid4 * w4 * (yVelocities[corner4] - prevYVelocities[corner4])) / d;
                float flipVelocity = velocity + weightedVelocityChanges;
                k2_y = PIC * picVelocity + (1 - PIC) * flipVelocity;
            }

            // Finally, derive the 3rd stage velocities
            float k3_x = 0;
            float k3_y = 0;
            float thirdXPosition = particleXPositions[i] + 0.75 * dt * k2_x;
            float thirdYPosition = particleYPositions[i] + 0.75 * dt * k2_y;

            // Interpolate the velocity field to get velocities at this position
            cellXCoordinate = floor(thirdXPosition / spacing);
            cellYCoordinate = floor((thirdYPosition - yShift) / spacing);
            // Extract the remainders (to compute corner weights w1, w2, w3, w4)
            deltaX = thirdXPosition - cellXCoordinate * spacing;
            deltaY = (thirdYPosition - yShift) - cellYCoordinate * spacing;
            // Now construct the corner weights
            w1 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing));
            w2 = (deltaX / spacing) * (1.0 - (deltaY / spacing));
            w3 = (deltaX / spacing) * (deltaY / spacing);
            w4 = (1.0 - (deltaX / spacing)) * (deltaY / spacing);
            // From here, compute the cell corners responsible for each weight
            corner1 = cellXCoordinate * gridHeight + cellYCoordinate;
            corner2 = fmin(cellXCoordinate + 1, gridLength - 2) * gridHeight + cellYCoordinate;
            corner3 = fmin(cellXCoordinate + 1, gridLength - 2) * gridHeight + fmin(cellYCoordinate + 1, gridHeight - 2);
            corner4 = cellXCoordinate * gridHeight + fmin(cellYCoordinate + 1, gridHeight - 2);
            // Proceed with transfer from the grid to particle
            offset = gridHeight;
            // Determine whether contributing corners on grid are valid or not
            valid1 = cellType[corner1] != EMPTY || cellType[corner1 - offset] != EMPTY ? 1.0 : 1.0;
            valid2 = cellType[corner2] != EMPTY || cellType[corner2 - offset] != EMPTY ? 1.0 : 1.0;
            valid3 = cellType[corner3] != EMPTY || cellType[corner3 - offset] != EMPTY ? 1.0 : 1.0;
            valid4 = cellType[corner4] != EMPTY || cellType[corner4 - offset] != EMPTY ? 1.0 : 1.0;
            d = w1 * valid1 + w2 * valid2 + w3 * valid3 + w4 * valid4;
            velocity = particleXVelocities[i];
            if (d > 0.0)    {
                float picVelocity = (valid1 * w1 * xVelocities[corner1] + valid2 * w2 * xVelocities[corner2] + valid3 * w3 * xVelocities[corner3] + valid4 * w4 * xVelocities[corner4]) / d;
                float weightedVelocityChanges = (valid1 * w1 * (xVelocities[corner1] - prevXVelocities[corner1]) + valid2 * w2 * (xVelocities[corner2] - prevXVelocities[corner2]) + valid3 * w3 * (xVelocities[corner3] - prevXVelocities[corner3]) + valid4 * w4 * (xVelocities[corner4] - prevXVelocities[corner4])) / d;
                float flipVelocity = velocity + weightedVelocityChanges;
                k3_x = PIC * picVelocity + (1 - PIC) * flipVelocity;
            }

            // Do the same for stage 3 y velocities
            cellXCoordinate = floor((thirdXPosition - xShift) / spacing);
            cellYCoordinate = floor(thirdYPosition / spacing);
            // Extract the remainders (to compute corner weights w1, w2, w3, w4)
            deltaX = (thirdXPosition - xShift) - cellXCoordinate * spacing;
            deltaY = thirdYPosition - cellYCoordinate * spacing;
            // Now construct the corner weights
            w1 = (1.0 - (deltaX / spacing)) * (1.0 - (deltaY / spacing));
            w2 = (deltaX / spacing) * (1.0 - (deltaY / spacing));
            w3 = (deltaX / spacing) * (deltaY / spacing);
            w4 = (1.0 - (deltaX / spacing)) * (deltaY / spacing);
            // From here, compute the cell corners responsible for each weight
            corner1 = cellXCoordinate * gridHeight + cellYCoordinate;
            corner2 = fmin(cellXCoordinate + 1, gridLength - 2) * gridHeight + cellYCoordinate;
            corner3 = fmin(cellXCoordinate + 1, gridLength - 2) * gridHeight + fmin(cellYCoordinate + 1, gridHeight - 2);
            corner4 = cellXCoordinate * gridHeight + fmin(cellYCoordinate + 1, gridHeight - 2);
            // Proceed with transfer from the grid to particle
            offset = 1;
            // Determine whether contributing corners on grid are valid or not
            valid1 = cellType[corner1] != EMPTY || cellType[corner1 - offset] != EMPTY ? 1.0 : 1.0;
            valid2 = cellType[corner2] != EMPTY || cellType[corner2 - offset] != EMPTY ? 1.0 : 1.0;
            valid3 = cellType[corner3] != EMPTY || cellType[corner3 - offset] != EMPTY ? 1.0 : 1.0;
            valid4 = cellType[corner4] != EMPTY || cellType[corner4 - offset] != EMPTY ? 1.0 : 1.0;
            d = w1 * valid1 + w2 * valid2 + w3 * valid3 + w4 * valid4;
            velocity = particleYVelocities[i];
            if (d > 0.0)    {
                float picVelocity = (valid1 * w1 * yVelocities[corner1] + valid2 * w2 * yVelocities[corner2] + valid3 * w3 * yVelocities[corner3] + valid4 * w4 * yVelocities[corner4]) / d;
                float weightedVelocityChanges = (valid1 * w1 * (yVelocities[corner1] - prevYVelocities[corner1]) + valid2 * w2 * (yVelocities[corner2] - prevYVelocities[corner2]) + valid3 * w3 * (yVelocities[corner3] - prevYVelocities[corner3]) + valid4 * w4 * (yVelocities[corner4] - prevYVelocities[corner4])) / d;
                float flipVelocity = velocity + weightedVelocityChanges;
                k3_y = PIC * picVelocity + (1 - PIC) * flipVelocity;
            }

            // Finally, update positions
            particleXPositions[i] += (2.0/9.0) * k1_x + (3.0/9.0) * k2_x + (5.0/9.0) * k3_x;
            particleYPositions[i] += (2.0/9.0) * k1_y + (3.0/9.0) * k2_y + (5.0/9.0) * k3_y;

            // Also update velocities
            particleYVelocities[i] += 9.81 * dt;
        }
    }
    else   {
        for (int i = 0; i < numParticles; i++)  {
            // Only velocity in the Y direction is integrated - due to gravity
            particleYVelocities[i] += dt * 9.81f;
            particleXPositions[i] += dt * particleXVelocities[i];
            particleYPositions[i] += dt * particleYVelocities[i];
        }
    }
}

void Fluid::extrapolateVelocities() {
    // extend the velocity field throughout the grid (not just where the particles are)
    // use an empty vector to store the 'wavefront' for BFS
    std::vector<std::tuple<int, int>> wX;
    std::vector<std::tuple<int, int>> wY;
    // Part 1: X extrapolation
    // loop over the entire grid (barring solid cells at grid edges)
    for (int i = 1; i < gridLength; i++)    {
        for (int j = 1; j < gridHeight; j++) {
            if (cellType[IX(i, j)] == SOLID) continue;
            bool neighborXMarker = (xMarker[IX(i + 1, j)] == 0 || xMarker[IX(i - 1, j)] == 0 || xMarker[IX(i, j + 1)] == 0 || xMarker[IX(i, j - 1)] == 0);
            if (xMarker[IX(i, j)] != 0 && neighborXMarker)  {
                xMarker[IX(i, j)] = 0;
                wX.emplace_back(i, j);
            }
        }
    }
    for (int t = 0; t < wX.size(); t++) {
        int i = std::get<0>(wX[t]);
        int j = std::get<1>(wX[t]);
        int numValidNeighbors = 0;
        int sum = 0;
        // check left neighbor
        if (cellType[IX(i - 1, j)] != SOLID) {
            if (xMarker[IX(i - 1, j)] < xMarker[IX(i, j)]) {
                numValidNeighbors++;
                sum += xVelocities[IX(i - 1, j)];
            }
            else if (xMarker[IX(i - 1, j)] == INT_MAX) {
                xMarker[IX(i - 1, j)] = xMarker[IX(i, j)] + 1;
                wX.push_back(std::tuple(i - 1, j));
            }
        }
        // check right neighbor
        if (cellType[IX(i + 1, j)] != SOLID) {
            if (xMarker[IX(i + 1, j)] < xMarker[IX(i, j)]) {
                numValidNeighbors++;
                sum += xVelocities[IX(i + 1, j)];
            }
            else if (xMarker[IX(i + 1, j)] == INT_MAX)  {
                xMarker[IX(i + 1, j)] = xMarker[IX(i, j)] + 1;
                wX.push_back(std::tuple(i + 1, j));
            }
        }
        // check bottom neighbor
        if (cellType[IX(i, j - 1)] != SOLID) {
            if (xMarker[IX(i, j - 1)] < xMarker[IX(i, j)]) {
                numValidNeighbors++;
                sum += xVelocities[IX(i, j - 1)];
            }
            else if (xMarker[IX(i, j - 1)] == INT_MAX)  {
                xMarker[IX(i, j - 1)] = xMarker[IX(i, j)] + 1;
                wX.push_back(std::tuple(i, j - 1));
            }
        }
        // check top neighbor
        if (cellType[IX(i, j + 1)] != SOLID) {
            if (xMarker[IX(i, j + 1)] < xMarker[IX(i, j)]) {
                numValidNeighbors++;
                sum += xVelocities[IX(i, j + 1)];
            }
            else if (xMarker[IX(i, j + 1)] == INT_MAX)  {
                xMarker[IX(i, j + 1)] = xMarker[IX(i, j)] + 1;
                wX.push_back(std::tuple(i, j + 1));
            }
        }
        if (numValidNeighbors > 0) xVelocities[IX(i, j)] = sum / numValidNeighbors;
    }

    // Part 2: Y extrapolation
    // loop over the entire grid (barring solid cells at grid edges)
    for (int i = 1; i < gridLength; i++)    {
        for (int j = 1; j < gridHeight; j++) {
            if (cellType[IX(i, j)] == SOLID) continue;
            bool neighborYMarker = (yMarker[IX(i + 1, j)] == 0 || yMarker[IX(i - 1, j)] == 0 || yMarker[IX(i, j + 1)] == 0 || yMarker[IX(i, j - 1)] == 0);
            if (yMarker[IX(i, j)] != 0 && neighborYMarker)  {
                yMarker[IX(i, j)] = 0;
                wY.emplace_back(i, j);
            }
        }
    }
    for (int t = 0; t < wY.size(); t++) {
        int i = std::get<0>(wY[t]);
        int j = std::get<1>(wY[t]);
        int numValidNeighbors = 0;
        int sum = 0;
        // check left neighbor
        if (cellType[IX(i - 1, j)] != SOLID) {
            if (yMarker[IX(i - 1, j)] < yMarker[IX(i, j)]) {
                numValidNeighbors++;
                sum += yVelocities[IX(i - 1, j)];
            }
            else if (yMarker[IX(i - 1, j)] == INT_MAX) {
                yMarker[IX(i - 1, j)] = yMarker[IX(i, j)] + 1;
                wY.push_back(std::tuple(i - 1, j));
            }
        }
        // check right neighbor
        if (cellType[IX(i + 1, j)] != SOLID) {
            if (yMarker[IX(i + 1, j)] < yMarker[IX(i, j)]) {
                numValidNeighbors++;
                sum += yVelocities[IX(i + 1, j)];
            }
            else if (yMarker[IX(i + 1, j)] == INT_MAX)  {
                yMarker[IX(i + 1, j)] = yMarker[IX(i, j)] + 1;
                wY.push_back(std::tuple(i + 1, j));
            }
        }
        // check bottom neighbor
        if (cellType[IX(i, j - 1)] != SOLID) {
            if (yMarker[IX(i, j - 1)] < yMarker[IX(i, j)]) {
                numValidNeighbors++;
                sum += yVelocities[IX(i, j - 1)];
            }
            else if (yMarker[IX(i, j - 1)] == INT_MAX)  {
                yMarker[IX(i, j - 1)] = yMarker[IX(i, j)] + 1;
                wY.push_back(std::tuple(i, j - 1));
            }
        }
        // check top neighbor
        if (cellType[IX(i, j + 1)] != SOLID) {
            if (yMarker[IX(i, j + 1)] < yMarker[IX(i, j)]) {
                numValidNeighbors++;
                sum += yVelocities[IX(i, j + 1)];
            }
            else if (yMarker[IX(i, j + 1)] == INT_MAX)  {
                yMarker[IX(i, j + 1)] = yMarker[IX(i, j)] + 1;
                wY.push_back(std::tuple(i, j + 1));
            }
        }
        if (numValidNeighbors > 0) yVelocities[IX(i, j)] = sum / numValidNeighbors;
    }
}

void Fluid::detectBoundaryCollisions() {
    // First, set a minimum and maximum limit for particles to be located at
    float minX = spacing + particleRadius;
    float minY = spacing + particleRadius;
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
