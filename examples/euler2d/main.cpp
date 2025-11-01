#include "../../src/hc.hpp"
#include <format>
#include <fstream>

using namespace hc;

extern void jump_ic(const GridGeo &geo, const IdealGasLaw &eos, GridArray &U);
extern void cylindrical_ic(const GridGeo &geo, const IdealGasLaw &eos, GridArray &U);
extern void four_piece_ic(const GridGeo &geo, const IdealGasLaw &eos, GridArray &U);

int main(int argc, char *argv[])
{
    // Simulation parameters
    const int Nx = 200;
    const int Ny = 200;
    const int ng = 1;
    const double Lx = 2.0;
    const double Ly = 2.0;
    const double t_start = 0.0;
    const double t_final = 0.25;
    const bool use_dim_split = true;
    const double CFL_1D = 0.7;
    const double CFL = use_dim_split ? CFL_1D : 0.5 * CFL_1D;

    // Equation of state and flux function
    IdealGasLaw eos(1.4);
    EulerFlux flux(&eos, 2);

    // Grid and geometry
    Grid grid(Nx, Ny, 1);
    GridGeo geo(grid, {0, 0, 0}, {Lx, Ly, 1});

    // State vectors
    GridArray U_new(grid, flux.n_comp(), {ng, ng, 0});
    GridArray U_old(grid, flux.n_comp(), {ng, ng, 0});

    // Initial condition
    cylindrical_ic(geo, eos, U_old);
    writeVTK("post/solution_0.vtk", geo, U_old, {"rho", "momx", "momy", "E"});

    // Compute initial timestep size
    double dt = compute_dt_estimate(U_old, geo, flux, CFL);

    // Numerical flux
    // ForceFlux nflux(&flux, &geo, &dt);
    HLLFlux nflux(&flux);
    // HLLCFlux nflux(&flux);

    // Maximum wavespeed
    double smax = -1;

    // Integrator
    auto erk = create_erk(U_old, create_rhs(geo, nflux, smax, true), ERKType::fe);

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

        if (use_dim_split)
        {
            const double smax_ = advance_split(geo, nflux, U_old, U_new, dt);
            if (smax_ > smax)
            {
                smax = smax_;
            }
        }
        else
        {
            for (int dir = 0; dir < flux.dim(); dir++)
            {
                U_old.fill_boundary(dir);
            }
            erk.advance(U_old, U_new, t, dt);
            GridArray::copy(U_new, U_old);
        }

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
        const double y = geo.cell_center(j, 1);
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
        const double x = geo.cell_center(i, 0)-1.0;
        const double y = geo.cell_center(j, 1)-1.0;
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

void four_piece_ic(const GridGeo &geo, const IdealGasLaw &eos, GridArray &U)
{
    const int dim = U.n_comp() - 2;
    const double gamma = eos.gamma();

    const double rho_high = 1.0;
    const double rho_medium = 0.4;
    const double rho_low = 0.125;

    const double p_high = 1.0;
    const double p_medium = 0.3;
    const double p_low = 0.1;

    const double x_lim = 0.75;
    const double y_lim = 0.75;

    grid_loop(U.grid(), [&](int i, int j, int k)
              {
        const double x = geo.cell_center(i, 0);
        const double y = geo.cell_center(j, 1);
        const double r = std::sqrt(x*x + y*y);
        const std::array<int, 3> idx = {i, j, k};
          if (x > x_lim)
        {
            if (y > y_lim)
            {
                U(idx, 0) = rho_high;
                U(idx, dim + 1) = p_high / (gamma - 1.0);
            }
            else
            {
                U(idx, 0) = rho_medium;
                U(idx, dim + 1) = p_medium / (gamma - 1.0);
            }
        }
        else
        {
            if (y > y_lim)
            {
                U(idx, 0) = rho_medium;
                U(idx, dim + 1) = p_medium / (gamma - 1.0);
            }
            else
            {
                U(idx, 0) = rho_low;
                U(idx, dim + 1) = p_low / (gamma - 1.0);
            }
        } });
}