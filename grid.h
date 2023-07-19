//
// Created by Anchit Mishra on 2023-07-19.
//

#ifndef FLUIDSIM_GRID_H
#define FLUIDSIM_GRID_H

#include <vector>

class grid {
public:
    // store any constant acceleration values
    float acceleration_x, acceleration_y;

    // store information about grid bounds
    int dim_x, dim_y; // number of grid cells - this is supposed to include walls on the outer layers
    float len_x, len_y; // actual length and width of the grid
    float spacing; // the length/width of one grid cell = len_x / dim_x
    float inv_spacing; // inverse of spacing, stored as an optimization

    // store information about fields stored on grid
    std::vector<float> u;
    std::vector<float> v;
    std::vector<float> u_prev;
    std::vector<float> v_prev;
    std::vector<char> cell_type;
    std::vector<float> phi; // for distance fields
    std::vector<double> pressure; // for PCG

    // primary functions to be performed on the grid
    void add_acceleration(); // -> to include accelerations in the velocity field
    void project(); // -> the main function to enforce incompressibility
    void extend_velocity_field(); // -> use fast sweeping to obtain valid velocity values outside fluid domain
    std::vector<float> sample_velocity_field(float x, float y) const;

    // utility functions to perform common math operations
    std::vector<float> compute_distance_offsets(const std::vector<float> &x, bool offset_y) const;
    std::vector<float> compute_bilerp_weights(const std::vector<float> &offsets) const;
    float bilerp(const std::vector<float> &pos, const std::vector<float> &f, const std::vector<float> &w) const;
    int IX(int x, int y) const; // -> return the cell coordinate in 1D indices
};


#endif //FLUIDSIM_GRID_H
