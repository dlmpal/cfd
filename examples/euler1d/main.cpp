#include "../../src/hc.hpp"
#include <format>
#include <fstream>

using namespace hc;

extern void jump_ic(const GridGeo &geo, const IdealGasLaw &eos, GridArray &U);
extern void save_timestamp(const std::string &filepath, const GridGeo &geo,
                           const GridArray &U);

int main(int argc, char *argv[])
{
    // Simulation parameters
    const int Nx = 200;
    const int ng = 1;
    const double Lx = 1.0;
    const double dx = Lx / (double)Nx;
    const double t_start = 0.0;
    const double t_final = 0.25;
    const double CFL = 0.8;

    // Equation of state and flux function
    IdealGasLaw eos(1.4);
    EulerFlux flux(&eos, 1);

    // Grid and geometry
    Grid grid(Nx, Nx, 1);
    GridGeo geo(grid, {0, 0, 0}, {Lx, 1, 1});

    // State and flux vectors
    GridArray U_new(grid, flux.n_comp(), {ng, 0, 0});
    GridArray U_old(grid, flux.n_comp(), {ng, 0, 0});
    GridArray Fx(grid.convert({1, 0, 0}), flux.n_comp());

    // Initial condition
    jump_ic(geo, eos, U_old);
    save_timestamp("post/solution_0", geo, U_old);
    writeVTK("post/solution_0.vtk", geo, U_old, {"rho", "momx", "E"});

    // Compute initial timestep size
    double dt = compute_dt_estimate(U_old, geo, flux, CFL);

    // Numerical flux
    // ForceFlux nflux(&flux, &dx, &dt);
    RusanovFlux nflux(&flux);

    // Maximum wavespeed
    double smax = -1;

    auto rhs_func = [&](const GridArray &S, GridArray &rhs, double time)
    {
        // Compute cell slopes
        GridArray slopes(S.grid(), S.n_comp(), {ng, 0, 0});
        compute_slopes(S, slopes, 0);
        slopes.fill_boundary(0);

        // Compute face fluxes
        const double smax_ = compute_fluxes(slopes, S, nflux, geo, Fx, dt, 0, false);
        if (smax_ > smax)
        {
            smax = smax_;
        }

        for (int n = 0; n < S.n_comp(); n++)
            grid_loop(S.grid(), [&](int i, int j, int k)
                      { rhs(i, j, k, n) = -1 / dx * (Fx(i + 1, j, k, n) - Fx(i, j, k, n)); });
    };
    auto erk = create_erk(U_old, rhs_func, ERKType::rk3);

    // Select times to save solution to file
    // const std::vector<double> timestamps = {t_final};
    // int timestamp_idx = 0;

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
        erk.advance(U_old, U_new, t, dt);
        GridArray::copy(U_new, U_old);

        // Save solution to file if required
        // if (std::abs(t - timestamps[timestamp_idx]) < 0.5 * dt)
        // {
        //     save_timestamp(std::format("post/solution_{}", timestamp_idx + 1), geo, U_new);
        //     timestamp_idx++;
        // }

        if (timestep % plot_int == 0)
        {
            writeVTK(std::format("post/solution_{}.vtk", timestep), geo, U_new, {"rho", "momx", "E"});
        }

        // Update timestep size
        dt = CFL * geo.dx() / smax;
    }

    return 0;
}

void jump_ic(const GridGeo &geo, const IdealGasLaw &eos, GridArray &U)
{
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

    const double dx = geo.dx();
    grid_loop(U.grid(), [&](int i, int j, int k)
              {
        const double x = i * dx;
        const std::array<int, 3> idx = {i, j, k};

        if (x < 0.5)
        {
            U(idx, 0) = rho_l;
            U(idx, 1) = rho_l * u_l;
            U(idx, 2) = E_l;
        }
        else
        {
            U(idx, 0) = rho_r;
            U(idx, 1) = rho_r * u_r;
            U(idx, 2) = E_r;
        } });
}

void save_timestamp(const std::string &filepath,
                    const GridGeo &geo,
                    const GridArray &U)
{
    std::ofstream file(filepath);
    const int Nx = U.grid().Nx();
    const double gamma = 1.4;
    for (int i = 0; i < Nx; i++)
    {
        const std::array<int, 3> idx = {i, 0, 0};
        const double rho = U(idx, 0);
        const double u = U(idx, 1) / rho;
        const double E = U(idx, 2);
        const double e = E / rho - 0.5 * u * u;
        const double p = (gamma - 1) * rho * e;
        file << i * geo.dx() << ", " << rho << ", " << u << ", " << p << ", " << e << "\n";
    }
}