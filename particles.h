//
// Created by Anchit Mishra on 2023-07-19.
//

#ifndef FLUIDSIM_PARTICLE_H
#define FLUIDSIM_PARTICLE_H

#include <vector>
#include "grid.h"

class particles {
    // store reference to grid
    grid& grid;

    // store positions and velocities
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> u;
    std::vector<float> v;

    // additional array to store grid normalization factors
    std::vector<float> acc;

    // store general information related to particles in the simulation
    int particleCount;

    // functions related to fluid simulation
    void spawn_particle(std::vector<float> pos, std::vector<float> vel); // -> for dynamically adding/deleting particles (keeping in mind the max allowed number of particles)
    void move_particles(float dt); // -> for moving particles given a valid velocity field
    void transfer_to_grid(); // -> for transferring particle velocities to grid nodes
    void transfer_from_grid(); // -> for copying corrected velocities from the grid
};


#endif //FLUIDSIM_PARTICLE_H
