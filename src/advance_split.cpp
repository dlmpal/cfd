#include "advance_split.hpp"
#include "flux_assembly.hpp"

namespace hc
{
    double advance_split(const GridGeo &geo,
                         const NumericalFlux &nflux,
                         GridArray &U_old,
                         GridArray &U_new,
                         double dt)
    {
        const FluxFunction &flux = *nflux.flux();
        const int dim = flux.dim();
        const int n_comp = flux.n_comp();

        const Grid &grid = U_old.grid();
        const auto [xlow, ylow, zlow] = grid.start();
        const auto [xhigh, yhigh, zhigh] = grid.end();

        double smax = -1;

        for (int dir = 0; dir < dim; dir++)
        {
            U_old.fill_boundary(dir);

            std::array<int, 3> n_ghost = {};
            n_ghost[dir] = U_old.n_ghost(dir);

            GridArray slopes(grid, n_comp, n_ghost);
            compute_slopes(U_old, slopes, dir);
            slopes.fill_boundary(dir);

            GridArray F(grid, n_comp, n_ghost);
            const double smax_ = compute_fluxes(slopes, U_old, nflux, geo, F, dt, dir, true);
            if (smax_ > smax)
            {
                smax = smax_;
            }

            for (int n = 0; n < n_comp; n++)
            {
                for (int k = zlow; k < zhigh; k++)
                {
                    for (int j = ylow; j < yhigh; j++)
                    {
                        for (int i = xlow; i < xhigh; i++)
                        {
                            std::array<int, 3> idx_left = {i, j, k};
                            std::array<int, 3> idx_right = {i, j, k};
                            idx_right[dir] += 1;
                            const double df = F(idx_right, n) - F(idx_left, n);
                            U_new(i, j, k, n) = U_old(i, j, k, n) - dt * geo.invdX(dir) * df;
                        }
                    }
                }
            }
            GridArray::copy(U_new, U_old);
        }

        return smax;
    }
}