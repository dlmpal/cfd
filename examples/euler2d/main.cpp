#include "../../src/hc.hpp"
#include <format>
#include <fstream>

using namespace hc;

extern void jump_ic(const GridGeo &geo, const IdealGasLaw &eos, GridArray &U);
extern void cylindrical_ic(const GridGeo &geo, const IdealGasLaw &eos, GridArray &U);

int main(int argc, char *argv[])
{
    // Simulation parameters
    const int Nx = 300;
    const int Ny = 300;
    const int ng = 1;
    const double Lx = 2.0;
    const double Ly = 2.0;
    const double dx = Lx / (double)Nx;
    const double dy = Ly / (double)Ny;
    const double t_start = 0.0;
    const double t_final = 0.25;
    const double CFL = 0.4;

    // Equation of state and flux function
    IdealGasLaw eos(1.4);
    EulerFlux flux(&eos, 2);

    // Grid and geometry
    Grid grid(Nx, Ny, 1);
    GridGeo geo(grid, {0, 0, 0}, {Lx, Ly, 1});

    // State and flux vectors
    GridArray U_new(grid, flux.n_comp(), {ng, ng, 0});
    GridArray U_old(grid, flux.n_comp(), {ng, ng, 0});
    std::array<GridArray, 2> F = {GridArray(grid.convert({1, 0, 0}), flux.n_comp()),
                                  GridArray(grid.convert({0, 1, 0}), flux.n_comp())};

    // Initial condition
    cylindrical_ic(geo, eos, U_old);
    writeVTK("post/solution_0.vtk", geo, U_old, {"rho", "momx", "momy", "E"});

    // Compute initial timestep size
    double dt = compute_dt_estimate(U_old, geo, flux, CFL);

    // Numerical flux
    ForceFlux nflux(&flux, &geo, &dt);

    // Maximum wavespeed
    double smax = -1;

    auto rhs_func = [&](const GridArray &S, GridArray &rhs, double time)
    {
        // Compute face fluxes (per direction)
        for (int dir = 0; dir < flux.dim(); dir++)
        {
            std::array<int, 3> n_ghost = {};
            n_ghost[dir] = ng;
            GridArray slopes(S.grid(), S.n_comp(), n_ghost);
            compute_slopes(S, slopes, dir);
            slopes.fill_boundary(dir);

            const double smax_ = compute_fluxes(slopes, S, nflux, geo, F[dir], dt, dir, false);
            if (smax_ > smax)
            {
                smax = smax_;
            }
        }

        rhs.set_all(0.0);

        for (int n = 0; n < S.n_comp(); n++)
        {
            grid_loop(S.grid(), [&](int i, int j, int k)
                      {
                        const double fx = F[0](i + 1, j, k, n) - F[0](i, j, k, n);
                        const double fy = F[1](i, j + 1, k, n) - F[1](i, j, k, n);
                        rhs(i, j, k, n) -= geo.invdX(0) * fx + geo.invdX(1) * fy; });
        }
    };
    auto erk = create_erk(U_old, rhs_func, ERKType::rk3);

    // Timestepping
    const int plot_int = 5;
    int timestep = 0;
    double t = t_start;
    while (t <= t_final)
    {
        // Increment time
        t += dt;
        timestep++;
        std::cout << std::format("Time: {}, Timestep: {}\n", t, timestep);

        U_old.fill_boundary(0);
        U_old.fill_boundary(1);

        erk.advance(U_old, U_new, t, dt);
        GridArray::copy(U_new, U_old);

        if (timestep % plot_int == 0)
        {
            writeVTK(std::format("post/solution_{}.vtk", timestep), geo, U_new, {"rho", "momx", "momy", "E"});
        }

        // Update timestep size
        dt = CFL * std::min({geo.dx(), geo.dy(), geo.dz()}) / smax;
    }

    return 0;
}

void jump_ic(const GridGeo &geo, const IdealGasLaw &eos, GridArray &U)
{
    const int dim = U.n_comp() - 2;
    const double dx = geo.dx();

    const double rho_l = 1.0;
    const double u_l = 0.0;
    const double p_l = 1.0;
    const double e_l = p_l / rho_l / (eos.gamma() - 1);
    const double E_l = rho_l * (e_l + 0.5 * u_l * u_l);

    const double rho_r = 0.125;
    const double u_r = 0.0;
    const double p_r = 0.1;
    const double e_r = p_r / rho_r / (eos.gamma() - 1);
    const double E_r = rho_r * (e_r + 0.5 * u_r * u_r);

    grid_loop(U.grid(), [&](int i, int j, int k)
              {
        const double y = j * geo.dy();
        const std::array<int, 3> idx = {i, j, k};

        if (y < 0.5)
        {
            U(idx, 0) = rho_l;
            U(idx, 2) = rho_l * u_l;
            U(idx, dim+1) = E_l;
        }
        else
        {
            U(idx, 0) = rho_r;
            U(idx, 2) = rho_r * u_r;
            U(idx, dim+1) = E_r;
        } });
}

void cylindrical_ic(const GridGeo &geo, const IdealGasLaw &eos, GridArray &U)
{
    const int dim = U.n_comp() - 2;
    const double dx = geo.dx();

    const double rho_l = 1.0;
    const double u_l = 0.0;
    const double v_l = 0.0;
    const double p_l = 1.0;
    const double e_l = p_l / rho_l / (eos.gamma() - 1);
    const double E_l = rho_l * e_l;

    const double rho_r = 0.125;
    const double u_r = 0.0;
    const double v_r = 0.0;
    const double p_r = 0.1;
    const double e_r = p_r / rho_r / (eos.gamma() - 1);
    const double E_r = rho_r * e_r;

    grid_loop(U.grid(), [&](int i, int j, int k)
              {
        const double x = i * geo.dx() - 1;
        const double y = j * geo.dy() - 1;
        const double r = std::sqrt(x*x + y*y);
        const std::array<int, 3> idx = {i, j, k};
        if (r < 0.4)
        {
            U(idx, 0) = rho_l;
            U(idx, 2) = rho_l * u_l;
            U(idx, dim+1) = E_l;
        }
        else
        {
            U(idx, 0) = rho_r;
            U(idx, 2) = rho_r * u_r;
            U(idx, dim+1) = E_r;
        } });
}