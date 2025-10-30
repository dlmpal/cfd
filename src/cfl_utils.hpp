#pragma once

#include "grid_array.hpp"
#include "flux_function.hpp"

namespace hc
{
    double compute_dt_estimate(const GridArray &U,
                               const GridGeo &geo,
                               const FluxFunction &flux,
                               double CFL);
}