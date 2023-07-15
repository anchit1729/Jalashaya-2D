//
// Created by Anchit Mishra on 2023-04-07.
//
#include "Fluid.h"
#define MAX_ITERS 200

void Fluid::projectPCG() {
    // Preconditioned Conjugate Gradient (Eigen implementation)
    SparseMatrix<double> A(numCells, numCells);
    VectorXd pressure(numCells);
    VectorXd divergence(numCells);
    prevXVelocities = xVelocities;
    prevYVelocities = yVelocities;
    A.setZero();
    divergence.setZero();
    pressure.setZero();
    A.reserve(VectorXd::Constant(numCells, 5)); // at most 5 non-zero elements per row in 2D
    // Step 1: Set up the divergence matrix
    float scale = 1.0 / spacing;
    for (int cellXCoordinate = 1; cellXCoordinate < gridLength-1; cellXCoordinate++)    {
        for (int cellYCoordinate = 1; cellYCoordinate < gridHeight-1; cellYCoordinate++)    {
            if (cellType[IX(cellXCoordinate, cellYCoordinate)] == FLUID)    {
                divergence[IX(cellXCoordinate, cellYCoordinate)] = -scale * (xVelocities[IX(cellXCoordinate+1, cellYCoordinate)] - xVelocities[IX(cellXCoordinate, cellYCoordinate)] + yVelocities[IX(cellXCoordinate, cellYCoordinate+1)] - yVelocities[IX(cellXCoordinate, cellYCoordinate)]);
            }
        }
    }
    for (int cellXCoordinate = 1; cellXCoordinate < gridLength-1; cellXCoordinate++)  {
        for (int cellYCoordinate = 1; cellYCoordinate < gridHeight-1; cellYCoordinate++)  {
            if (cellType[cellXCoordinate * gridHeight + cellYCoordinate] == FLUID)  {
                // compute the negative divergence for all fluid cells
                int cellCoordinate = cellXCoordinate * gridHeight + cellYCoordinate;
                int rightCellCoordinate = (cellXCoordinate + 1) * gridHeight + cellYCoordinate;
                int topCellCoordinate = cellXCoordinate * gridHeight + cellYCoordinate + 1;
                int leftCellCoordinate = (cellXCoordinate - 1) * gridHeight + cellYCoordinate;
                int bottomCellCoordinate = cellXCoordinate * gridHeight + cellYCoordinate - 1;
                // When accounting for solids, we assume solid wall velocity = 0
                if (cellType[leftCellCoordinate] == SOLID)    {
                    divergence[cellCoordinate] -= scale * (xVelocities[cellCoordinate]);
                }
                if (cellType[rightCellCoordinate] == SOLID) {
                    divergence[cellCoordinate] += scale * (xVelocities[rightCellCoordinate]);
                }
                if (cellType[bottomCellCoordinate] == SOLID)    {
                    divergence[cellCoordinate] -= scale * (yVelocities[cellCoordinate]);
                }
                if (cellType[topCellCoordinate] == SOLID)   {
                    divergence[cellCoordinate] += scale * (yVelocities[topCellCoordinate]);
                }
            }
        }
    }
    // Step 2: Set up the A matrix entries
    // NOTE: We exploit symmetry, thus using only positive direction variables
    float dt = TIMESTEP / SUBSTEPS;
    scale = dt / (density * spacing * spacing);
    for (int cellXCoordinate = 1; cellXCoordinate < gridLength-1; cellXCoordinate++)  {
        for (int cellYCoordinate = 1; cellYCoordinate < gridHeight-1; cellYCoordinate++)  {
            if (cellType[cellXCoordinate * gridHeight + cellYCoordinate] == FLUID)  {
                int leftCellCoordinate = (cellXCoordinate - 1) * gridHeight + cellYCoordinate;
                int rightCellCoordinate = (cellXCoordinate + 1) * gridHeight + cellYCoordinate;
                int bottomCellCoordinate = cellXCoordinate * gridHeight + cellYCoordinate - 1;
                int topCellCoordinate = cellXCoordinate * gridHeight + cellYCoordinate + 1;
                int cellCoordinate = cellXCoordinate * gridHeight + cellYCoordinate;
                // handle negative x neighbor
                if (cellType[leftCellCoordinate] == FLUID)  {
                    A.coeffRef(cellCoordinate, cellCoordinate) += scale;
                }
                // handle positive x neighbor
                if (cellType[rightCellCoordinate] == FLUID) {
                    A.coeffRef(cellCoordinate, cellCoordinate) += scale;
                    A.coeffRef(cellCoordinate, rightCellCoordinate) = -scale;
                    A.coeffRef(rightCellCoordinate, cellCoordinate) = -scale;
                }
                else if (cellType[rightCellCoordinate] == EMPTY)    {
                    A.coeffRef(cellCoordinate, cellCoordinate) += scale;
                }
                // handle negative y neighbor
                if (cellType[bottomCellCoordinate] == FLUID)    {
                    A.coeffRef(cellCoordinate, cellCoordinate) += scale;
                }
                // handle positive y neighbor
                if (cellType[topCellCoordinate] == FLUID)   {
                    A.coeffRef(cellCoordinate, cellCoordinate) += scale;
                    A.coeffRef(cellCoordinate, topCellCoordinate) = -scale;
                    A.coeffRef(topCellCoordinate, cellCoordinate) = -scale;
                }
                else if (cellType[topCellCoordinate] == EMPTY)  {
                    A.coeffRef(cellCoordinate, cellCoordinate) += scale;
                }
            }
        }
    }
    // Step 3: Run the PCG algorithm with Eigen
    ConjugateGradient<SparseMatrix<double>, Eigen::UpLoType::Lower|Eigen::UpLoType::Upper> cg;
    cg.setMaxIterations(250);
    cg.setTolerance(1e-6);
    cg.compute(A);
    pressure = cg.solve(divergence);
    std::cout << "Iteration count: " << cg.iterations() << std::endl;
    std::cout << "Estimated error: " << cg.error() << std::endl;
    // Step 4: Use pressure estimates to compute velocity updates
    scale = dt / (density * spacing);
    for (int cellXCoordinate = 1; cellXCoordinate < gridLength-1; cellXCoordinate++)  {
        for (int cellYCoordinate = 1; cellYCoordinate < gridHeight-1; cellYCoordinate++)  {
            int leftCellCoordinate = (cellXCoordinate - 1) * gridHeight + cellYCoordinate;
            int bottomCellCoordinate = cellXCoordinate * gridHeight + cellYCoordinate - 1;
            int cellCoordinate = cellXCoordinate * gridHeight + cellYCoordinate;
            if (cellType[leftCellCoordinate] == FLUID || cellType[cellCoordinate] == FLUID) {
                if (cellType[leftCellCoordinate] == SOLID || cellType[cellCoordinate] == SOLID) {
                    xVelocities[cellCoordinate] = 0;
                } else  {
                    xVelocities[cellCoordinate] -= scale * (pressure[cellCoordinate] - pressure[leftCellCoordinate]);
                }
            } else  {
                // mark unknown
                xVelocities[cellCoordinate] = 0;
            }
            if (cellType[bottomCellCoordinate] == FLUID || cellType[cellCoordinate] == FLUID)   {
                if (cellType[bottomCellCoordinate] == SOLID || cellType[cellCoordinate] == SOLID)   {
                    yVelocities[cellCoordinate] = 0;
                } else  {
                    yVelocities[cellCoordinate] -= scale * (pressure[cellCoordinate] - pressure(bottomCellCoordinate));
                }
            } else  {
                // mark unknown
                yVelocities[cellCoordinate] = 0;
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
                    // Compute multiplication factor (inspired by Ten Minute Physics' explanation to account for solid cells)
                    float leftMultiplicationFactor = cellType[leftCellCoordinate] != SOLID ? 1.0 : 0.0;
                    float rightMultiplicationFactor = cellType[rightCellCoordinate] != SOLID ? 1.0 : 0.0;
                    float topMultiplicationFactor = cellType[topCellCoordinate] != SOLID ? 1.0 : 0.0;
                    float bottomMultiplicationFactor = cellType[bottomCellCoordinate] != SOLID ? 1.0 : 0.0;
                    float normalizationFactor = leftMultiplicationFactor + rightMultiplicationFactor + topMultiplicationFactor + bottomMultiplicationFactor;
                    leftMultiplicationFactor /= normalizationFactor;
                    rightMultiplicationFactor /= normalizationFactor;
                    topMultiplicationFactor /= normalizationFactor;
                    bottomMultiplicationFactor /= normalizationFactor;
                    // Finally, compute divergence and modify velocities considering all multiplication factors
                    divergence += (rightXVelocity - leftXVelocity + topYVelocity - bottomYVelocity);
                    divergence *= overRelaxation;
                    // Account for the particle density by subtracting the compression factor
                    float stiffnessCoefficient = 1;
                    if (initialParticleDensity > 0)    {
                        float compression = cellParticleDensity[cellCoordinate] - initialParticleDensity;
                        if (compression > 0) divergence -= (stiffnessCoefficient * compression);// / spacing;
                    }
                    xVelocities[cellCoordinate] += divergence * leftMultiplicationFactor;
                    xVelocities[rightCellCoordinate] -= divergence * rightMultiplicationFactor;
                    yVelocities[cellCoordinate] += divergence * bottomMultiplicationFactor;
                    yVelocities[topCellCoordinate] -= divergence * topMultiplicationFactor;
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
