#include "grid_array.hpp"

namespace hc
{
    //=============================================================================
    static std::array<int, 4> compute_stides(int Nx, int Ny, int Nz)
    {
        std::array<int, 4> strides;
        strides[0] = 1.0;
        strides[1] = strides[0] * Nx;
        strides[2] = strides[1] * Ny;
        strides[3] = strides[2] * Nz;
        return strides;
    }
    //=============================================================================
    GridArray::GridArray(const Grid &grid, int n_comp, const std::array<int, 3> &n_ghost)
        : grid_(grid),
          n_comp_(n_comp),
          n_ghost_(n_ghost),
          stride_(compute_stides(grid.Nx() + 2 * n_ghost_[0],
                                 grid.Ny() + 2 * n_ghost_[1],
                                 grid.Nz() + 2 * n_ghost_[2])),
          values_(stride_[3] * n_comp, 0.0)
    {
    }
    //=============================================================================
    GridArray::GridArray(const Grid &grid, int n_comp, int n_ghost)
        : GridArray(grid, n_comp, {n_ghost, n_ghost, n_ghost})
    {
    }
    //=============================================================================
    Grid GridArray::grid() const
    {
        return grid_;
    }
    //=============================================================================
    int GridArray::n_comp() const
    {
        return n_comp_;
    }
    //=============================================================================
    int GridArray::n_ghost(int dir) const
    {
        return n_ghost_[dir];
    }
    //=============================================================================
    int GridArray::flattened_index(int i, int j, int k, int n) const
    {
        i += n_ghost_[0];
        j += n_ghost_[1];
        k += n_ghost_[2];
        const int idx = i + stride_[1] * j + stride_[2] * k + stride_[3] * n;
        if (idx < 0 || idx >= stride_[3] * n_comp_)
        {
            std::cerr << "Invalid index: " << idx << "\n";
            exit(EXIT_FAILURE);
        }
        return idx;
    }
    //=============================================================================
    double GridArray::operator()(int i, int j, int k, int n) const
    {
        return values_[flattened_index(i, j, k, n)];
    }
    //=============================================================================
    double &GridArray::operator()(int i, int j, int k, int n)
    {
        return values_[flattened_index(i, j, k, n)];
    }
    //=============================================================================
    double GridArray::operator()(const std::array<int, 3> &idx, int n) const
    {
        return values_[flattened_index(idx[0], idx[1], idx[2], n)];
    }
    //=============================================================================
    double &GridArray::operator()(const std::array<int, 3> &idx, int n)
    {
        return values_[flattened_index(idx[0], idx[1], idx[2], n)];
    }
    //=============================================================================
    void GridArray::set_all(double value)
    {
        std::fill(values_.begin(), values_.end(), value);
    }
    //=============================================================================
    void GridArray::fill_boundary(int dir)
    {
        const int Nx = grid_.Nx();
        const int Ny = grid_.Ny();
        const int Nz = grid_.Nz();

        for (int n = 0; n < n_comp_; n++)
        {
            grid_loop(grid_, [&](int i, int j, int k)
                      {
                          if (dir == 0)
                          {
                              (*this)(-1, j, k, n) = (*this)(0, j, k, n);
                              (*this)(Nx, j, k, n) = (*this)(Nx - 1, j, k, n);
                          }
                          else if (dir == 1)
                          {
                              (*this)(i, -1, k, n) = (*this)(i, 0, k, n);
                              (*this)(i, Ny, k, n) = (*this)(i, Ny - 1, k, n);
                          }
                          else if (dir == 2)
                          {
                              (*this)(i, j, -1, n) = (*this)(i, j, 0, n);
                              (*this)(i, j, Nz, n) = (*this)(i, j, Nz-1, n);
                          } });
        }
    }
    //=============================================================================
    void GridArray::copy(const GridArray &src, GridArray &dest)
    {
        std::copy(src.values_.cbegin(), src.values_.cend(), dest.values_.begin());
    }
    //=============================================================================
    void GridArray::axpy(double a, const GridArray &x, GridArray &y)
    {
        for (int n = 0; n < x.n_comp_; n++)
        {
            grid_loop(x.grid_, [&](int i, int j, int k)
                      { y(i, j, k, n) += a * x(i, j, k, n); });
        }
    }
}