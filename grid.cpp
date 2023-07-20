//
// Created by Anchit Mishra on 2023-07-19.
//

#include "grid.h"

int grid::IX(int x, int y) const   {
    return x * dim_y + y;
}

std::vector<float> grid::compute_bilerp_weights(const std::vector<float> &offsets) const {
    // note: order assumes starting at top left of a cell, working clockwise
    // note: y-axis increases in the downward direction
    std::vector<float> w = {0, 0, 0, 0};
    w[0] = (1 - (offsets[0] * inv_spacing)) * (1 - (offsets[1] * inv_spacing));
    w[1] = (offsets[0] * inv_spacing) * (1 - (offsets[1] * inv_spacing));
    w[2] = (offsets[0] * inv_spacing) * (offsets[1] * inv_spacing);
    w[3] = (1 - (offsets[0] * inv_spacing)) * (offsets[1] * inv_spacing);
    return w;
}

std::vector<float> grid::compute_distance_offsets(const std::vector<float> &x, const bool offset_y = true) const {
    float x0 = x[0], y0 = x[1];
    float x_shift = offset_y ? 0.0f : static_cast<float>(0.5 * spacing);
    float y_shift = offset_y ? static_cast<float>(0.5 * spacing) : 0.0f;
    x0 -= x_shift; std::clamp(x0, spacing, static_cast<float>((dim_x - 2)) * spacing); // bound within domain
    y0 -= y_shift; std::clamp(y0, spacing, static_cast<float>((dim_y - 2)) * spacing); // bound within domain
    int cell_x = floor(x0 * inv_spacing); std::clamp<int>(cell_x, 1, dim_x - 2); // bound within domain
    int cell_y = floor(y0 * inv_spacing); std::clamp<int>(cell_y, 1, dim_y - 2); // bound within domain
    float dx = x0 - static_cast<float>(cell_x) * spacing, dy = y0 - static_cast<float>(cell_y) * spacing;
    std::vector<float> offsets = {dx, dy};
    return offsets;
}

std::vector<float> grid::sample_velocity_field(float x, float y) const {
    std::vector<float> pos = {x, y};
    std::vector<float> offsets;
    offsets = compute_distance_offsets(pos, true);
    std::vector<float> w = compute_bilerp_weights(offsets);
    std::vector<float> vel;
    pos = {x, static_cast<float>(y - 0.5 * spacing)};
    float res1 = grid::bilerp(pos, u, w);

    offsets = compute_distance_offsets(pos, false);
    w = compute_bilerp_weights(offsets);
    pos = {static_cast<float>(x - 0.5 * spacing), y};
    float res2 = grid::bilerp(pos, v, w);

    return std::vector<float>{res1, res2};
}

float grid::bilerp(const std::vector<float> &pos, const std::vector<float> &f, const std::vector<float> &w) const {
    float val = 0;
    int x = static_cast<int>(floorf(pos[0] * inv_spacing));
    int y = static_cast<int>(floorf(pos[1] * inv_spacing));
    int c1 = IX(x, y);
    int c2 = IX(x + 1, y);
    int c3 = IX(x + 1, y + 1);
    int c4 = IX(x, y + 1);
    std::vector<int> c = {c1, c2, c3, c4};
    for (int i = 0; i < 4; i++) {
        val += f[c[i]] * w[i];
    }
    return val;
}
