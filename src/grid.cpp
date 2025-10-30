#include "grid.hpp"

namespace hc
{
    //=============================================================================
    Grid::Grid(int Nx, int Ny, int Nz, const std::array<int, 3> nodality)
        : Nx_(Nx),
          Ny_(Ny),
          Nz_(Nz),
          nodality_(nodality)
    {
        /// @todo Check nodality
    }
    //=============================================================================
    int Grid::Nx() const
    {
        return Nx_ + nodality_[0];
    }
    //=============================================================================
    int Grid::Ny() const
    {
        return Ny_ + nodality_[1];
    }
    //=============================================================================
    int Grid::Nz() const
    {
        return Nz_ + nodality_[2];
    }
    //=============================================================================
    Grid Grid::convert(const std::array<int, 3> nodality) const
    {
        return Grid(Nx_, Ny_, Nz_, nodality);
    }
    //=============================================================================
    std::array<int, 3> Grid::start() const
    {
        return {0, 0, 0};
    }
    //=============================================================================
    std::array<int, 3> Grid::end() const
    {
        return {Nx_ + nodality_[0],
                Ny_ + nodality_[1],
                Nz_ + nodality_[2]};
    }
    //=============================================================================
    GridGeo::GridGeo(const Grid &grid,
                     const std::array<double, 3> &start,
                     const std::array<double, 3> &end)
        : grid_(grid),
          start_(start),
          end_(end),
          dX_({(end[0] - start[0]) / (double)grid.Nx(),
               (end[1] - start[1]) / (double)grid.Ny(),
               (end[2] - start[2]) / (double)grid.Nz()}),
          invdX_({1 / dX_[0],
                  1 / dX_[1],
                  1 / dX_[2]})
    {
    }
    //=============================================================================
    double GridGeo::start(int dir) const
    {
        return start_[dir];
    }
    //=============================================================================
    double GridGeo::end(int dir) const
    {
        return end_[dir];
    }
    //=============================================================================
    double GridGeo::dX(int dir) const
    {
        return dX_[dir];
    }
    //=============================================================================
    double GridGeo::invdX(int dir) const
    {
        return invdX_[dir];
    }
    //=============================================================================
    double GridGeo::dx() const
    {
        return dX_[0];
    }
    //=============================================================================
    double GridGeo::dy() const
    {
        return dX_[1];
    }
    //=============================================================================
    double GridGeo::dz() const
    {
        return dX_[2];
    }
}