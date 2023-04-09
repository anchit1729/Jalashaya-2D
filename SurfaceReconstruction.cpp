//
// Created by Anchit Mishra on 2023-04-08.
//
#include "Fluid.h"

#define SWEEP_ITERS 2

// Implementing fast sweeping to construct a signed distance field
float Fluid::distance(std::tuple<int, int> x, std::tuple<float, float> p)   {
    int x1 = std::get<0>(x);
    int x2 = std::get<1>(x);
    float p1 = std::get<0>(p);
    float p2 = std::get<1>(p);
    return sqrt(pow(p1 - x1, 2) + pow(p2 - x2, 2));
}

void Fluid::computeSignedDistance() {
    std::fill(pointDistance.begin(), pointDistance.end(), INT_MAX);
    std::fill(closestPointIndices.begin(), closestPointIndices.end(), -1);
    // step 1: initialize the arrays near the geometry
    // loop over the particles
    for (int i = 0; i < numParticles; i++)  {
        // locate the point's grid cell coordinates
        float x = particleXPositions[i];
        float y = particleYPositions[i];
        int lowX = floor(x);
        int lowY = floor(y);
        // check distance
        if (distance(std::tuple(lowX, lowY), std::tuple(x, y)) < pointDistance[IX(floor(x / spacing), floor(y / spacing))]) {
            pointDistance[IX(floor(x / spacing), floor(y / spacing))] = distance(std::tuple(lowX, lowY), std::tuple(x, y));
            closestPointIndices[IX(floor(x / spacing), floor(y / spacing))] = i;
        }
    }

    // step 2: propagate distances
    // order 1: ascending, ascending
    for (int i = 0; i < gridLength; i++)    {
        for (int j = 0; j < gridHeight; j++)    {
            int currCoord = IX(i, j);
            int leftX = fmax(i - 1, 0);
            int leftY = j;
            int leftCoord = IX(leftX, leftY);
            int rightX = fmin(i + 1, gridLength - 1);
            int rightY = j;
            int rightCoord = IX(rightX, rightY);
            int bottomX = i;
            int bottomY = fmax(j - 1, 0);
            int bottomCoord = IX(bottomX, bottomY);
            int topX = i;
            int topY = fmin(j + 1, gridHeight - 1);
            int topCoord = IX(topX, topY);
            if (leftCoord != currCoord && closestPointIndices[leftCoord] != -1) {
                // consider this point
                if (distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[leftCoord]], particleYPositions[closestPointIndices[leftCoord]])) < pointDistance[currCoord])  {
                    pointDistance[currCoord] = distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[leftCoord]], particleYPositions[closestPointIndices[leftCoord]]));
                    closestPointIndices[currCoord] = closestPointIndices[leftCoord];
                }
            }
            if (rightCoord != currCoord && closestPointIndices[rightCoord] != -1) {
                // consider this point
                if (distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[rightCoord]], particleYPositions[closestPointIndices[rightCoord]])) < pointDistance[currCoord])  {
                    pointDistance[currCoord] = distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[rightCoord]], particleYPositions[closestPointIndices[rightCoord]]));
                    closestPointIndices[currCoord] = closestPointIndices[rightCoord];
                }
            }
            if (topCoord != currCoord && closestPointIndices[topCoord] != -1) {
                // consider this point
                if (distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[topCoord]], particleYPositions[closestPointIndices[topCoord]])) < pointDistance[currCoord])  {
                    pointDistance[currCoord] = distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[topCoord]], particleYPositions[closestPointIndices[topCoord]]));
                    closestPointIndices[currCoord] = closestPointIndices[topCoord];
                }
            }
            if (bottomCoord != currCoord && closestPointIndices[bottomCoord] != -1) {
                // consider this point
                if (distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[bottomCoord]], particleYPositions[closestPointIndices[bottomCoord]])) < pointDistance[currCoord])  {
                    pointDistance[currCoord] = distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[bottomCoord]], particleYPositions[closestPointIndices[bottomCoord]]));
                    closestPointIndices[currCoord] = closestPointIndices[bottomCoord];
                }
            }
        }
    }
    // order 2: ascending, descending
    for (int i = 0; i < gridLength; i++)    {
        for (int j = gridHeight - 1; j >= 0; j--)   {
            int currCoord = IX(i, j);
            int leftX = fmax(i - 1, 0);
            int leftY = j;
            int leftCoord = IX(leftX, leftY);
            int rightX = fmin(i + 1, gridLength - 1);
            int rightY = j;
            int rightCoord = IX(rightX, rightY);
            int bottomX = i;
            int bottomY = fmax(j - 1, 0);
            int bottomCoord = IX(bottomX, bottomY);
            int topX = i;
            int topY = fmin(j + 1, gridHeight - 1);
            int topCoord = IX(topX, topY);
            if (leftCoord != currCoord && closestPointIndices[leftCoord] != -1) {
                // consider this point
                if (distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[leftCoord]], particleYPositions[closestPointIndices[leftCoord]])) < pointDistance[currCoord])  {
                    pointDistance[currCoord] = distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[leftCoord]], particleYPositions[closestPointIndices[leftCoord]]));
                    closestPointIndices[currCoord] = closestPointIndices[leftCoord];
                }
            }
            if (rightCoord != currCoord && closestPointIndices[rightCoord] != -1) {
                // consider this point
                if (distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[rightCoord]], particleYPositions[closestPointIndices[rightCoord]])) < pointDistance[currCoord])  {
                    pointDistance[currCoord] = distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[rightCoord]], particleYPositions[closestPointIndices[rightCoord]]));
                    closestPointIndices[currCoord] = closestPointIndices[rightCoord];
                }
            }
            if (topCoord != currCoord && closestPointIndices[topCoord] != -1) {
                // consider this point
                if (distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[topCoord]], particleYPositions[closestPointIndices[topCoord]])) < pointDistance[currCoord])  {
                    pointDistance[currCoord] = distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[topCoord]], particleYPositions[closestPointIndices[topCoord]]));
                    closestPointIndices[currCoord] = closestPointIndices[topCoord];
                }
            }
            if (bottomCoord != currCoord && closestPointIndices[bottomCoord] != -1) {
                // consider this point
                if (distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[bottomCoord]], particleYPositions[closestPointIndices[bottomCoord]])) < pointDistance[currCoord])  {
                    pointDistance[currCoord] = distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[bottomCoord]], particleYPositions[closestPointIndices[bottomCoord]]));
                    closestPointIndices[currCoord] = closestPointIndices[bottomCoord];
                }
            }
        }
    }
    // order 3: descending, ascending
    for (int i = gridLength - 1; i >= 0; i--)    {
        for (int j = 0; j < gridHeight; j++)    {
            int currCoord = IX(i, j);
            int leftX = fmax(i - 1, 0);
            int leftY = j;
            int leftCoord = IX(leftX, leftY);
            int rightX = fmin(i + 1, gridLength - 1);
            int rightY = j;
            int rightCoord = IX(rightX, rightY);
            int bottomX = i;
            int bottomY = fmax(j - 1, 0);
            int bottomCoord = IX(bottomX, bottomY);
            int topX = i;
            int topY = fmin(j + 1, gridHeight - 1);
            int topCoord = IX(topX, topY);
            if (leftCoord != currCoord && closestPointIndices[leftCoord] != -1) {
                // consider this point
                if (distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[leftCoord]], particleYPositions[closestPointIndices[leftCoord]])) < pointDistance[currCoord])  {
                    pointDistance[currCoord] = distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[leftCoord]], particleYPositions[closestPointIndices[leftCoord]]));
                    closestPointIndices[currCoord] = closestPointIndices[leftCoord];
                }
            }
            if (rightCoord != currCoord && closestPointIndices[rightCoord] != -1) {
                // consider this point
                if (distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[rightCoord]], particleYPositions[closestPointIndices[rightCoord]])) < pointDistance[currCoord])  {
                    pointDistance[currCoord] = distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[rightCoord]], particleYPositions[closestPointIndices[rightCoord]]));
                    closestPointIndices[currCoord] = closestPointIndices[rightCoord];
                }
            }
            if (topCoord != currCoord && closestPointIndices[topCoord] != -1) {
                // consider this point
                if (distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[topCoord]], particleYPositions[closestPointIndices[topCoord]])) < pointDistance[currCoord])  {
                    pointDistance[currCoord] = distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[topCoord]], particleYPositions[closestPointIndices[topCoord]]));
                    closestPointIndices[currCoord] = closestPointIndices[topCoord];
                }
            }
            if (bottomCoord != currCoord && closestPointIndices[bottomCoord] != -1) {
                // consider this point
                if (distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[bottomCoord]], particleYPositions[closestPointIndices[bottomCoord]])) < pointDistance[currCoord])  {
                    pointDistance[currCoord] = distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[bottomCoord]], particleYPositions[closestPointIndices[bottomCoord]]));
                    closestPointIndices[currCoord] = closestPointIndices[bottomCoord];
                }
            }
        }
    }
    // order 4: descending, descending
    for (int i = gridLength - 1; i >= 0; i--)   {
        for (int j = gridHeight - 1; j >= 0; j--)   {
            int currCoord = IX(i, j);
            int leftX = fmax(i - 1, 0);
            int leftY = j;
            int leftCoord = IX(leftX, leftY);
            int rightX = fmin(i + 1, gridLength - 1);
            int rightY = j;
            int rightCoord = IX(rightX, rightY);
            int bottomX = i;
            int bottomY = fmax(j - 1, 0);
            int bottomCoord = IX(bottomX, bottomY);
            int topX = i;
            int topY = fmin(j + 1, gridHeight - 1);
            int topCoord = IX(topX, topY);
            if (leftCoord != currCoord && closestPointIndices[leftCoord] != -1) {
                // consider this point
                if (distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[leftCoord]], particleYPositions[closestPointIndices[leftCoord]])) < pointDistance[currCoord])  {
                    pointDistance[currCoord] = distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[leftCoord]], particleYPositions[closestPointIndices[leftCoord]]));
                    closestPointIndices[currCoord] = closestPointIndices[leftCoord];
                }
            }
            if (rightCoord != currCoord && closestPointIndices[rightCoord] != -1) {
                // consider this point
                if (distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[rightCoord]], particleYPositions[closestPointIndices[rightCoord]])) < pointDistance[currCoord])  {
                    pointDistance[currCoord] = distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[rightCoord]], particleYPositions[closestPointIndices[rightCoord]]));
                    closestPointIndices[currCoord] = closestPointIndices[rightCoord];
                }
            }
            if (topCoord != currCoord && closestPointIndices[topCoord] != -1) {
                // consider this point
                if (distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[topCoord]], particleYPositions[closestPointIndices[topCoord]])) < pointDistance[currCoord])  {
                    pointDistance[currCoord] = distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[topCoord]], particleYPositions[closestPointIndices[topCoord]]));
                    closestPointIndices[currCoord] = closestPointIndices[topCoord];
                }
            }
            if (bottomCoord != currCoord && closestPointIndices[bottomCoord] != -1) {
                // consider this point
                if (distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[bottomCoord]], particleYPositions[closestPointIndices[bottomCoord]])) < pointDistance[currCoord])  {
                    pointDistance[currCoord] = distance(std::tuple(i, j), std::tuple(particleXPositions[closestPointIndices[bottomCoord]], particleYPositions[closestPointIndices[bottomCoord]]));
                    closestPointIndices[currCoord] = closestPointIndices[bottomCoord];
                }
            }
        }
    }
}