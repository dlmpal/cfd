#pragma once

#include "grid_array.hpp"
#include "numerical_flux.hpp"

namespace hc
{
    /// @brief
    /// @param U
    /// @param slopes
    /// @param dir
    /// @param w
    void compute_slopes(const GridArray &U,
                        GridArray &slopes,
                        int dir, double w = 0.0);
    /// @brief
    /// @param slopes
    /// @param U
    /// @param nflux
    /// @param geo
    /// @param F
    /// @param dt
    /// @param use_halfstep
    /// @return
    double compute_fluxes(const GridArray &slopes,
                          const GridArray &U,
                          const NumericalFlux &nflux,
                          const GridGeo &geo,
                          GridArray &F,
                          double dt, int dir,
                          bool use_halfstep = true);
}