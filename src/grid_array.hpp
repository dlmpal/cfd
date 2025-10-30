#pragma once

#include "grid.hpp"
#include <vector>
#include <iostream>

namespace hc
{
    class GridArray
    {
    public:
        GridArray(const Grid &grid, int n_comp, const std::array<int, 3> &n_ghost);
        GridArray(const Grid &grid, int n_comp = 1, int n_ghost = 0);

        Grid grid() const;
        int n_comp() const;
        int n_ghost(int dir = 0) const;

        double operator()(int i, int j, int k, int n) const;
        double &operator()(int i, int j, int k, int n);
        double operator()(const std::array<int, 3> &idx, int n) const;
        double &operator()(const std::array<int, 3> &idx, int n);

        void set_all(double value);

        /// @todo Implement
        void fill_periodic(int dir);

        /// @todo Allow for arbitrary number of ghosts
        void fill_boundary(int dir);

        static void copy(const GridArray &src, GridArray &dest);
        static void axpy(double a, const GridArray &x, GridArray &y);

    private:
        int flattened_index(int i, int j, int k, int n) const;

    private:
        const Grid grid_;

        const int n_comp_;

        const std::array<int, 3> n_ghost_;

        const std::array<int, 4> stride_;

        std::vector<double> values_;
    };
}