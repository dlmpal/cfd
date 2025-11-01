#pragma once

#include <array>
#include <concepts>

namespace hc
{
    class Grid
    {
    public:
        Grid(int Nx, int Ny, int Nz, const std::array<int, 3> nodality = {});

        int Nx() const;
        int Ny() const;
        int Nz() const;

        Grid convert(const std::array<int, 3> nodality) const;

        std::array<int, 3> start() const;
        std::array<int, 3> end() const;

    private:
        const int Nx_;
        const int Ny_;
        const int Nz_;
        const std::array<int, 3> nodality_;
    };

    class GridGeo
    {
    public:
        GridGeo(const Grid &grid,
                const std::array<double, 3> &start,
                const std::array<double, 3> &end);

        double start(int dir) const;
        double end(int dir) const;

        double dX(int dir) const;
        double invdX(int dir) const;

        double dx() const;
        double dy() const;
        double dz() const;

        double cell_center(int idx, int dir) const;

    private:
        const Grid grid_;
        const std::array<double, 3> start_;
        const std::array<double, 3> end_;
        const std::array<double, 3> dX_;
        const std::array<double, 3> invdX_;
    };

    template <typename F>
    concept GridFunc = requires(F f, int i, int j, int k) {
        { f(i, j, k) } -> std::same_as<void>;
    };

    void grid_loop(const Grid &grid, GridFunc auto &&func)
    {
        const auto start = grid.start();
        const auto end = grid.end();

        for (int k = start[2]; k < end[2]; k++)
        {
            for (int j = start[1]; j < end[1]; j++)
            {
                for (int i = start[0]; i < end[0]; i++)
                {
                    func(i, j, k);
                }
            }
        }
    }
}