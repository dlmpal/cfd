#include "erk.hpp"
#include "flux_assembly.hpp"

namespace hc
{
    //=============================================================================
    ERKIntegrator::ERKIntegrator(const GridArray &state, RHSFunction rhs,
                                 int n_stages, std::vector<double> &&nodes,
                                 std::vector<double> &&weights,
                                 std::vector<std::vector<double>> &&coeffs)
        : rhs_(rhs),
          n_stages_(n_stages),
          nodes_(std::move(nodes)),
          weights_(std::move(weights)),
          coeffs_(std::move(coeffs))
    {
        for (int i = 0; i < n_stages_; i++)
        {
            stages_.emplace_back(state.grid(), state.n_comp());
        }
    }
    //=============================================================================
    void ERKIntegrator::advance(const GridArray &S_old, GridArray &S_new, double time, double dt)
    {
        // Evaluate the RHS for each stage
        for (int i = 0; i < n_stages_; i++)
        {
            // Evaluate the intermediate state
            GridArray::copy(S_old, S_new);
            for (int j = 0; j < i; j++)
            {
                double aij = coeffs_[i][j];
                GridArray::axpy(dt * aij, stages_[j], S_new);
            }
            rhs_(S_new, stages_[i], time + dt * nodes_[i], dt);
        }

        // Evaluate the next state
        GridArray::copy(S_old, S_new);
        for (int i = 0; i < n_stages_; i++)
        {
            GridArray::axpy(dt * weights_[i], stages_[i], S_new);
        }
    }
    //=============================================================================
    ERKIntegrator create_fe(const GridArray &state, ERKIntegrator::RHSFunction rhs)
    {
        const int n_stages = 1;
        std::vector<double> nodes = {0.};
        std::vector<double> weights = {1.0};
        std::vector<std::vector<double>> coeffs = {{0.0}};

        return ERKIntegrator(state, rhs, n_stages,
                             std::move(nodes),
                             std::move(weights),
                             std::move(coeffs));
    }
    //=============================================================================
    ERKIntegrator create_rk2(const GridArray &state, ERKIntegrator::RHSFunction rhs)
    {
        const int n_stages = 2;
        const double a = 2. / 3.;
        std::vector<double> nodes = {0., a};
        std::vector<double> weights = {1. - 0.5 / a, 0.5 / a};
        std::vector<std::vector<double>> coeffs(n_stages);
        coeffs[0] = {0., 0.};
        coeffs[1] = {a, 0.};

        return ERKIntegrator(state, rhs, n_stages,
                             std::move(nodes),
                             std::move(weights),
                             std::move(coeffs));
    }
    //=============================================================================
    ERKIntegrator create_rk3(const GridArray &state, ERKIntegrator::RHSFunction rhs)
    {
        const int n_stages = 3;
        std::vector<double> nodes = {0, 1., 0.5};
        std::vector<double> weights = {1. / 6., 1. / 6., 2. / 3.};
        std::vector<std::vector<double>> coeffs(n_stages);
        coeffs[0] = {0., 0., 0.};
        coeffs[1] = {1., 0., 0.};
        coeffs[2] = {0.25, 0.25, 0.};

        return ERKIntegrator(state, rhs, n_stages,
                             std::move(nodes),
                             std::move(weights),
                             std::move(coeffs));
    }
    //=============================================================================
    ERKIntegrator create_rk4(const GridArray &state, ERKIntegrator::RHSFunction rhs)
    {
        const int n_stages = 4;
        std::vector<double> nodes = {0., 1. / 2., 1. / 2., 1.};
        std::vector<double> weights = {1. / 6., 1. / 3., 1. / 3., 1. / 6.};
        std::vector<std::vector<double>> coeffs(n_stages);
        coeffs[0] = {0., 0., 0., 0.};
        coeffs[1] = {1. / 2., 0., 0., 0.};
        coeffs[2] = {0., 1. / 2., 0., 0.};
        coeffs[3] = {0., 0., 1., 0.};

        return ERKIntegrator(state, rhs, n_stages,
                             std::move(nodes),
                             std::move(weights),
                             std::move(coeffs));
    }
    //=============================================================================
    ERKIntegrator create_erk(const GridArray &state, ERKIntegrator::RHSFunction rhs, ERKType type)
    {
        switch (type)
        {
        case ERKType::fe:
            return create_fe(state, rhs);
        case ERKType::rk2:
            return create_rk2(state, rhs);
        case ERKType::rk3:
            return create_rk3(state, rhs);
        case ERKType::rk4:
            return create_rk4(state, rhs);
        default:
            // SFEM_ERROR("Invalid ERK type\n");
            return ERKIntegrator(state, rhs, 0, {}, {}, {});
        }
    }
    //=============================================================================
    ERKIntegrator::RHSFunction create_rhs(const GridGeo &geo,
                                          const NumericalFlux &nflux,
                                          double &smax,
                                          bool use_halfstep)
    {
        return [&, use_halfstep](const GridArray &U, GridArray &rhs, double time, double dt)
        {
            rhs.set_all(0.0);

            const FluxFunction &flux = *nflux.flux();
            const int dim = flux.dim();
            const int n_comp = flux.n_comp();

            const Grid &grid = U.grid();
            const auto [xlow, ylow, zlow] = grid.start();
            const auto [xhigh, yhigh, zhigh] = grid.end();

            // Compute face fluxes (per direction)
            for (int dir = 0; dir < flux.dim(); dir++)
            {
                std::array<int, 3> n_ghost = {};
                n_ghost[dir] = U.n_ghost(dir);

                GridArray slopes(grid, n_comp, n_ghost);
                compute_slopes(U, slopes, dir);
                slopes.fill_boundary(dir);

                GridArray F(grid, n_comp, n_ghost);
                const double smax_ = compute_fluxes(slopes, U, nflux, geo, F, dt, dir, use_halfstep);
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
                                rhs(i, j, k, n) -= geo.invdX(dir) * df;
                            }
                        }
                    }
                }
            }
        };
    }
}