#pragma once

#include "numerical_flux.hpp"
#include "erk.hpp"

namespace hc
{
    double advance_split(const GridGeo &geo,
                         const NumericalFlux &nflux,
                         GridArray &U_old,
                         GridArray &U_new,
                         double dt);
}