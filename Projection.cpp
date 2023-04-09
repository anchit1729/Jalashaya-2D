//
// Created by Anchit Mishra on 2023-04-07.
//
#include "Fluid.h"
#define MAX_ITERS 200

void Fluid::projectPCG() {
    // Preconditioned Conjugate Gradient (vectorized, using Eigen)
    // Compute a mapping for fluid indices
    std::vector<int> fluidIndices;
    fluidIndices.resize(numCells, -1);
    int numFluidCells = 0;
    for (int i = 0; i < gridLength; i++)    {
        for (int j = 0; j < gridHeight; j++)    {
            if (cellType[IX(i, j)] == FLUID)    {
                fluidIndices[IX(i, j)] = numFluidCells;
                numFluidCells++;
            }
        }
    }
    // first, we need to prepare the sparse matrix A
    A.resize(numFluidCells, numFluidCells);
    pressure.resize(numFluidCells);
    divergence.resize(numFluidCells);
    A.setZero();
    pressure.setZero();
    divergence.setZero();
    A.reserve(VectorXd::Constant(numCells, 5)); // reserve room for 5 non-zero elements per row
    double scale = (TIMESTEP / SUBSTEPS) / (density * spacing * spacing);
    for (int i = 0; i < gridLength; i++)    {
        for (int j = 0; j < gridHeight; j++)    {
            if (cellType[IX(i, j)] == FLUID)    {
                // only populate entries for fluid cells
                // check neighboring cells to populate
                if (cellType[IX(i - 1, j)] == FLUID)    {
                    // left neighbor
                    A.coeffRef(fluidIndices[IX(i, j)], fluidIndices[IX(i - 1, j)]) = -scale;
                    A.coeffRef(fluidIndices[IX(i, j)], fluidIndices[IX(i, j)]) += scale;
                }
                else if (cellType[IX(i - 1, j)] == EMPTY)   {
                    A.coeffRef(fluidIndices[IX(i, j)], fluidIndices[IX(i, j)]) += scale;
                }
                if (cellType[IX(i + 1, j)] == FLUID)    {
                    // right neighbor
                    A.coeffRef(fluidIndices[IX(i, j)], fluidIndices[IX(i + 1, j)]) = -scale;
                    A.coeffRef(fluidIndices[IX(i, j)], fluidIndices[IX(i, j)]) += scale;
                }
                else if (cellType[IX(i + 1, j)] == EMPTY)   {
                    A.coeffRef(fluidIndices[IX(i, j)], fluidIndices[IX(i, j)]) += scale;
                }
                if (cellType[IX(i, j - 1)] == FLUID)    {
                    // bottom neighbor
                    A.coeffRef(fluidIndices[IX(i, j)], fluidIndices[IX(i, j - 1)]) = -scale;
                    A.coeffRef(fluidIndices[IX(i, j)], fluidIndices[IX(i, j)]) += scale;
                }
                else if (cellType[IX(i, j - 1)] == EMPTY)   {
                    A.coeffRef(fluidIndices[IX(i, j)], fluidIndices[IX(i, j)]) += scale;
                }
                if (cellType[IX(i, j + 1)] == FLUID)    {
                    // top neighbor
                    A.coeffRef(fluidIndices[IX(i, j)], fluidIndices[IX(i, j + 1)]) = -scale;
                    A.coeffRef(fluidIndices[IX(i, j)], fluidIndices[IX(i, j)]) += scale;
                }
                else if (cellType[IX(i, j + 1)] == EMPTY)   {
                    A.coeffRef(fluidIndices[IX(i, j)], fluidIndices[IX(i, j)]) += scale;
                }
            }
        }
    }

    // next, compute divergence for all fluid cells
    for (int i = 1; i < gridLength - 1; i++)    {
        for (int j = 1; j < gridHeight - 1; j++)    {
            if (cellType[IX(i, j)] == FLUID)    {
                // only compute for fluid cells
                double horizontal = (xVelocities[IX(i + 1, j)] - xVelocities[IX(i, j)]);
                double vertical = (yVelocities[IX(i, j + 1)] - yVelocities[IX(i, j)]);
                divergence.coeffRef(fluidIndices[IX(i, j)]) = -1 * (horizontal + vertical) / spacing;
            }
        }
    }
    // modify it to account for solid-fluid boundaries
    for (int i = 1; i < gridLength - 1; i++)    {
        for (int j = 1; j < gridHeight - 1; j++)    {
            if (cellType[IX(i, j)] == FLUID)    {
                // only for fluid cells
                if (cellType[IX(i - 1, j)] == SOLID)    {
                    divergence.coeffRef(fluidIndices[IX(i, j)]) -= (1 / spacing) * (xVelocities[IX(i, j)] - containerWallXVelocity);
                }
                if (cellType[IX(i + 1, j)] == SOLID)    {
                    divergence.coeffRef(fluidIndices[IX(i, j)]) += (1 / spacing) * (xVelocities[IX(i, j)] - containerWallXVelocity);
                }
                if (cellType[IX(i, j - 1)] == SOLID)    {
                    divergence.coeffRef(fluidIndices[IX(i, j)]) -= (1 / spacing) * (yVelocities[IX(i, j)] - containerWallYVelocity);
                }
                if (cellType[IX(i, j + 1)] == SOLID)    {
                    divergence.coeffRef(fluidIndices[IX(i, j)]) += (1 / spacing) * (yVelocities[IX(i, j)] - containerWallYVelocity);
                }
            }
        }
    }
    // system prepared, now precondition and solve!
    ConjugateGradient<SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cg;
    cg.setMaxIterations(MAX_ITERS);
    cg.setTolerance(1e-6);
    cg.compute(A);
    pressure = cg.solve(divergence);
    std::cout << "Iterations: " << cg.iterations() << "\n";
    std::cout << "Error: " << cg.error() << "\n";
    // update the velocity field based on the solved pressure values
    scale = (TIMESTEP / SUBSTEPS) / (density * spacing);
    for (int i = 1; i < gridLength - 1; i++)    {
        for (int j = 1; j < gridHeight - 1; j++)    {
            // first, the x component
            if (cellType[IX(i - 1, j)] == FLUID || cellType[IX(i, j)] == FLUID) {
                if (cellType[IX(i - 1, j)] == SOLID || cellType[IX(i, j)] == SOLID) {
                    xVelocities[IX(i, j)] = 0; xMarker[IX(i, j)] = 0; // solid velocities
                } else  {
                    float pressureDifference = 0;
                    int x = fluidIndices[IX(i, j)];
                    int y = fluidIndices[IX(i - 1, j)];
                    if (x != -1 && y != -1) pressureDifference = (x - y);
                    else if (x == -1 && y != -1) pressureDifference = -y;
                    else if (x != -1 && y == -1) pressureDifference = x;
                    else pressureDifference = 0;
                    xVelocities[IX(i, j)] -= scale * pressureDifference;
                    xMarker[IX(i, j)] = 0;
                }
            } else  {
                xVelocities[IX(i, j)] = 0;
                xMarker[IX(i, j)] = INT_MAX;
            }
            // next, the y component
            if (cellType[IX(i, j - 1)] == FLUID || cellType[IX(i, j)] == FLUID) {
                if (cellType[IX(i, j - 1)] == SOLID || cellType[IX(i, j)] == SOLID) {
                    yVelocities[IX(i, j)] = 0; yMarker[IX(i, j)] = 0; // solid velocities
                } else  {
                    float pressureDifference = 0;
                    int x = fluidIndices[IX(i, j)];
                    int y = fluidIndices[IX(i, j - 1)];
                    if (x != -1 && y != -1) pressureDifference = (x - y);
                    else if (x == -1 && y != -1) pressureDifference = -y;
                    else if (x != -1 && y == -1) pressureDifference = x;
                    else pressureDifference = 0;
                    yVelocities[IX(i, j)] -= scale * pressureDifference;
                    yMarker[IX(i, j)] = 0;
                }
            } else  {
                yVelocities[IX(i, j)] = 0;
                yMarker[IX(i, j)] = INT_MAX;
            }
        }
    }
}

void Fluid::projectGS() {
    // Gauss-Seidel solver (without vectorization)
    prevXVelocities = xVelocities;
    prevYVelocities = yVelocities;
    int numIterations = 100;
    float overRelaxation = 1.9;
    for (int i = 0; i < numIterations; i++) {
        // Note that we go from 1 to gridDimension - 2 to make sure there are always valid velocity components to pull from
        for (int cellXCoordinate = 1; cellXCoordinate < gridLength - 1; cellXCoordinate++)    {
            for (int cellYCoordinate = 1; cellYCoordinate < gridHeight - 1; cellYCoordinate++)    {
                if (cellType[cellXCoordinate * gridHeight + cellYCoordinate] == FLUID)  {
                    // Define top, bottom, left and right cell coordinates
                    int cellCoordinate = cellXCoordinate * gridHeight + cellYCoordinate;
                    int leftCellCoordinate = (cellXCoordinate - 1) * gridHeight + cellYCoordinate;
                    int rightCellCoordinate = (cellXCoordinate + 1) * gridHeight + cellYCoordinate;
                    int topCellCoordinate = cellXCoordinate * gridHeight + cellYCoordinate + 1;
                    int bottomCellCoordinate = cellXCoordinate * gridHeight + cellYCoordinate - 1;
                    // Extract velocity components
                    float leftXVelocity = xVelocities[cellCoordinate];
                    float rightXVelocity = xVelocities[rightCellCoordinate];
                    float bottomYVelocity = yVelocities[cellCoordinate];
                    float topYVelocity = yVelocities[topCellCoordinate];
                    // Initialize divergence
                    float divergence = 0.0;
                    // Account for the particle density by subtracting the compression factor
                    float stiffnessCoefficient = 1;
                    if (initialParticleDensity > 0)    {
                        float compression = cellParticleDensity[cellCoordinate] - initialParticleDensity;
                        if (compression > 0) divergence -= stiffnessCoefficient * compression;
                    }
                    // Compute multiplication factor (inspired by Ten Minute Physics' explanation to account for solid cells)
                    float leftMultiplicationFactor = cellType[leftCellCoordinate] == FLUID || cellType[leftCellCoordinate] == EMPTY ? 1.0 : 0.0;
                    float rightMultiplicationFactor = cellType[rightCellCoordinate] == FLUID || cellType[rightCellCoordinate] == EMPTY ? 1.0 : 0.0;
                    float topMultiplicationFactor = cellType[topCellCoordinate] == FLUID || cellType[topCellCoordinate] == EMPTY ? 1.0 : 0.0;
                    float bottomMultiplicationFactor = cellType[bottomCellCoordinate] == FLUID || cellType[bottomCellCoordinate] == EMPTY ? 1.0 : 0.0;
                    float normalizationFactor = leftMultiplicationFactor + rightMultiplicationFactor + topMultiplicationFactor + bottomMultiplicationFactor;
                    leftMultiplicationFactor /= normalizationFactor;
                    rightMultiplicationFactor /= normalizationFactor;
                    topMultiplicationFactor /= normalizationFactor;
                    bottomMultiplicationFactor /= normalizationFactor;
                    // Finally, compute divergence and modify velocities considering all multiplication factors
                    divergence += rightXVelocity - leftXVelocity + topYVelocity - bottomYVelocity;
                    xVelocities[cellCoordinate] += divergence * leftMultiplicationFactor * overRelaxation;
                    xVelocities[rightCellCoordinate] -= divergence * rightMultiplicationFactor * overRelaxation;
                    yVelocities[cellCoordinate] += divergence * bottomMultiplicationFactor * overRelaxation;
                    yVelocities[topCellCoordinate] -= divergence * topMultiplicationFactor * overRelaxation;
                }
            }
        }
    }
}

void Fluid::computeCellDensities() {
    // First, reset the current cellParticleDensity grid
    std::fill(cellParticleDensity.begin(), cellParticleDensity.end(), 0.0);
    // Since density is stored at cell centers, we need to shift both coordinates by half the grid spacing
    float xShift = 0.5f * spacing;
    float yShift = 0.5f * spacing;
    // Now, iterate over all particles
    for (int i = 0; i < numParticles; i++)  {
        // Same idea as velocity transfer - find out which cell a particle belongs to and update the particle density for that cell
        float x = Utils::clamp(particleXPositions[i], spacing, (gridLength - 1) * spacing);
        float y = Utils::clamp(particleYPositions[i], spacing, (gridHeight - 1) * spacing);
        int cellXCoordinate = floor((x - xShift) / spacing);
        int cellYCoordinate = floor((y - yShift) / spacing);
        float deltaX = (x - xShift) - (float)cellXCoordinate * spacing;
        float deltaY = (y - yShift) - (float)cellYCoordinate * spacing;
        // Compute the corner weights
        float w1 = (1 - deltaX / spacing) * (1 - deltaY / spacing);
        float w2 = (deltaX / spacing) * (1 - deltaY / spacing);
        float w3 = (deltaX / spacing) * (deltaY / spacing);
        float w4 = (1 - deltaX / spacing) * (deltaY / spacing);
        // Compute the corner coordinates
        int corner1 = cellXCoordinate * gridHeight + cellYCoordinate;
        int corner2 = fmin((float)cellXCoordinate + 1, (float)gridLength - 2) * gridHeight + cellYCoordinate;
        int corner3 = fmin(cellXCoordinate + 1, gridLength - 2) * gridHeight + fmin(cellYCoordinate + 1, gridHeight - 2);
        int corner4 = cellXCoordinate * gridHeight + fmin(cellYCoordinate + 1, gridHeight - 2);
        // Accumulate weight values at cellParticleDensity
        cellParticleDensity[corner1] += w1;
        cellParticleDensity[corner2] += w2;
        cellParticleDensity[corner3] += w3;
        cellParticleDensity[corner4] += w4;
    }

    // There is one more case to solve - when the fluid has just been initialized
    if (initialParticleDensity == 0.0) {
        float totalDensityAccumulator = 0;
        int fluidCellCount = 0;
        for (int i = 0; i < numCells; i++)  {
            if (cellType[i] == FLUID) {
                fluidCellCount++;
                totalDensityAccumulator += cellParticleDensity[i];
            }
        }
        if (fluidCellCount > 0) initialParticleDensity = totalDensityAccumulator / fluidCellCount;
    }
}
