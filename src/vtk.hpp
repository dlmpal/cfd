#pragma once

#include "grid_array.hpp"

namespace hc
{
    void writeVTK(const std::string &filename,
                  const GridGeo &geo,
                  const GridArray &arr,
                  const std::vector<std::string> &components);
}