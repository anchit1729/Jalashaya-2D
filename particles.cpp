//
// Created by Anchit Mishra on 2023-07-19.
//

#include "particles.h"

/**
 * Utility function for dynamically spawning new particles in a simulation with specified position and velocity.
 *
 * @param x - x-component of position
 * @param y - y-component of position
 * @param u - x-component of velocity
 * @param v - y-component of velocity
 */
void particles::spawn_particle(std::vector<float> pos, std::vector<float> vel) {
    this->x.push_back(pos[0]);
    this->y.push_back(pos[1]);
    this->u.push_back(vel[0]);
    this->v.push_back(vel[1]);
    particleCount++;
}

/**
 * Advection function - handles RK2 advection.
 *
 * @param dt - time-step size
 */
void particles::move_particles(float dt) {
    // define walls for clamping
    auto min_x = grid.spacing * 1.001, min_y = grid.spacing * 1.001;
    auto max_x = grid.len_x - grid.spacing * 1.001, max_y = grid.len_y - grid.spacing * 1.001;
    // RK2 advection
    for (int i = 0; i < particleCount; i++) {
        // step 1: compute intermediate step
        auto x1 = x[i] + 0.5 * dt * u[i];
        auto y1 = y[i] + 0.5 * dt * v[i];
        // step 2: clamp intermediate values
        std::clamp(x1, min_x, max_x);
        std::clamp(y1, min_y, max_y);
        // step 3: compute final values
        std::vector<float> v2 = grid.sample_velocity_field(static_cast<float>(x1), static_cast<float>(y1));
        float x2 = x[i] + dt * v2[0];
        float y2 = y[i] + dt * v2[1];
        // step 4: clamp final values
        std::clamp<float>(x2, static_cast<float>(min_x), static_cast<float>(max_x));
        std::clamp<float>(y2, static_cast<float>(min_y), static_cast<float>(max_y));
        x[i] = x2;
        y[i] = y2;
    }
}

/**
 * Splatting function - copies velocities from particle onto the grid.
 */
void particles::transfer_to_grid() {
    // step 1: reset all relevant variables
    grid.u_prev = grid.u;
    grid.v_prev = grid.v;
    std::fill(grid.u.begin(), grid.u.end(), 0.0);
    std::fill(grid.v.begin(), grid.v.end(), 0.0);
    std::fill(acc.begin(), acc.end(), 0.0);

    // step 2: transfer x components
    // step 3: transfer y components
    // step 4: mark cells that contain fluid
}

void particles::transfer_from_grid() {

}
