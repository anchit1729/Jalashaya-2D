//
// Created by Anchit Mishra on 2023-04-07.
//
#include "Fluid.h"

void Fluid::transferVelocitiesToGrid() {
    // Transfer velocities of particles to their nearby grid side faces
    // First, store a copy of previous velocity field values (for FLIP)
    prevXVelocities = xVelocities;
    prevYVelocities = yVelocities;
    // Now, begin with transferring all x velocity components
    for (int i = 0; i < numParticles; i++)  {
        // determine the 
    }
}