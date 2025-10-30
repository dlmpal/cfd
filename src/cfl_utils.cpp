#include "cfl_utils.hpp"
#include <algorithm>

namespace hc
{
    //=============================================================================
    double compute_dt_estimate(const GridArray &U,
                               const GridGeo &geo,
                               const FluxFunction &flux,
                               double CFL)
    {
        const Grid &grid = U.grid();
        const auto [xlow, ylow, zlow] = grid.start();
        const auto [xhigh, yhigh, zhigh] = grid.end();
        const int n_comp = U.n_comp();
        const int dim = flux.dim();

        std::vector<double> u(n_comp);
        std::vector<double> f(n_comp);

        // Maximum wavespeed
        double smax = -1;

        // Loop over all cells
        for (int k = zlow; k < zhigh; k++)
        {
            for (int j = ylow; j < yhigh; j++)
            {
                for (int i = xlow; i < xhigh; i++)
                {
                    // Get the state values for the cell
                    for (int n = 0; n < n_comp; n++)
                    {
                        u[n] = U(i, j, k, n);
                    }

                    // Evaluate the cell fluxes (per direction)
                    // Update the maximum wavespeed if needed
                    for (int dir = 0; dir < dim; dir++)
                    {
                        const double s = flux(u, f, dir);
                        if (s > smax)
                        {
                            smax = s;
                        }
                    }
                }
            }
        }

        return CFL * std::min({geo.dx(), geo.dy(), geo.dz()}) / smax;
    }
}