#include "flux_assembly.hpp"

namespace hc
{
    //=============================================================================
    static double minbee(double r)
    {
        if (r <= 0)
        {
            return 0;
        }
        else if (r < 1)
        {
            return r;
        }
        else
        {
            return std::min(1.0, 2 / (1 + r));
        }
    }
    //=============================================================================
    void compute_slopes(const GridArray &U,
                        GridArray &slopes,
                        int dir, double w)
    {
        const int n_comp = U.n_comp();
        const Grid &grid = U.grid();
        const auto [xlow, ylow, zlow] = grid.start();
        const auto [xhigh, yhigh, zhigh] = grid.end();

        for (int k = zlow; k < zhigh; k++)
        {
            for (int j = ylow; j < yhigh; j++)
            {
                for (int i = xlow; i < xhigh; i++)
                {
                    const std::array<int, 3> idx = {i, j, k};

                    std::array<int, 3> idx_l = {i, j, k};
                    idx_l[dir] -= 1;

                    std::array<int, 3> idx_r = {i, j, k};
                    idx_r[dir] += 1;

                    for (int n = 0; n < n_comp; n++)
                    {
                        const double delta_l = U(idx, n) - U(idx_l, n);
                        const double delta_r = U(idx_r, n) - U(idx, n);
                        const double delta = 0.5 * ((1 + w) * delta_l + (1 - w) * delta_r);
                        const double r = delta_l / delta_r;
                        slopes(idx, n) = minbee(r) * delta;
                    }
                }
            }
        }
    }
    //=============================================================================
    static void evolve_local(const std::array<int, 3> idx,
                             const GridArray &slopes,
                             const GridArray &U,
                             const FluxFunction &flux,
                             std::span<double> u,
                             const GridGeo &geo,
                             double dt,
                             int dir, double sign)
    {
        std::vector<double> ul(flux.n_comp());
        std::vector<double> ur(flux.n_comp());
        for (int i = 0; i < flux.n_comp(); i++)
        {
            ul[i] = U(idx, i) - 0.5 * slopes(idx, i);
            ur[i] = U(idx, i) + 0.5 * slopes(idx, i);
        }

        std::vector<double> fl(flux.n_comp());
        flux(ul, fl, dir);

        std::vector<double> fr(flux.n_comp());
        flux(ur, fr, dir);

        for (int i = 0; i < flux.n_comp(); i++)
        {
            const double df = fr[i] - fl[i];
            const double u_ = U(idx, i) + sign * 0.5 * slopes(idx, i);
            u[i] = u_ - 0.5 * dt / geo.dX(dir) * df;
        }
    }
    //=============================================================================
    double compute_fluxes(const GridArray &slopes,
                          const GridArray &U,
                          const NumericalFlux &nflux,
                          const GridGeo &geo,
                          GridArray &F,
                          double dt, int dir,
                          bool use_halfstep)
    {
        const FluxFunction &flux = *nflux.flux();
        const int dim = flux.dim();
        const int n_comp = flux.n_comp();

        std::vector<double> ul(n_comp);
        std::vector<double> ur(n_comp);
        std::vector<double> f(n_comp);

        // Maximum wave speed
        double smax = -1;

        std::array<int, 3> nodality = {};
        nodality[dir] = 1;

        const Grid grid = U.grid().convert(nodality);
        const auto [xlow, ylow, zlow] = grid.start();
        const auto [xhigh, yhigh, zhigh] = grid.end();

        for (int k = zlow; k < zhigh; k++)
        {
            for (int j = ylow; j < yhigh; j++)
            {
                for (int i = xlow; i < xhigh; i++)
                {
                    std::array<int, 3> idx_r = {i, j, k};
                    std::array<int, 3> idx_l = {i, j, k};
                    idx_l[dir] -= 1;

                    if (use_halfstep)
                    {
                        evolve_local(idx_l, slopes, U, flux, ul, geo, dt, dir, 1);
                        evolve_local(idx_r, slopes, U, flux, ur, geo, dt, dir, -1);
                    }
                    else
                    {
                        for (int n = 0; n < flux.n_comp(); n++)
                        {
                            ul[n] = U(idx_l, n) + 0.5 * slopes(idx_l, n);
                            ur[n] = U(idx_r, n) - 0.5 * slopes(idx_r, n);
                        }
                    }

                    const double s = nflux.compute_flux(ul, ur, f, dir);

                    for (int n = 0; n < flux.n_comp(); n++)
                    {
                        F(i, j, k, n) = f[n];
                    }

                    if (s > smax)
                    {
                        smax = s;
                    }
                }
            }
        }

        return smax;
    }
}