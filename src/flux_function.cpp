#include "flux_function.hpp"
#include "equation_of_state.hpp"

namespace hc
{
    //=============================================================================
    FluxFunction::FluxFunction(int dim, int n_comp)
        : dim_(dim),
          n_comp_(n_comp)
    {
    }
    //=============================================================================
    int FluxFunction::dim() const
    {
        return dim_;
    }
    //=============================================================================
    int FluxFunction::n_comp() const
    {
        return n_comp_;
    }
    //=============================================================================
    AdvectionFlux::AdvectionFlux(double speed)
        : FluxFunction(1, 1),
          speed_(speed)
    {
    }
    //=============================================================================
    double AdvectionFlux::speed() const
    {
        return speed_;
    }
    //=============================================================================
    double AdvectionFlux::operator()(std::span<const double> u, std::span<double> f, int dir) const
    {
        f[0] = speed_ * u[0];
        return speed_;
    }
    //=============================================================================
    BurgersFlux::BurgersFlux()
        : FluxFunction(1, 1)
    {
    }
    //=============================================================================
    double BurgersFlux::operator()(std::span<const double> u, std::span<double> f, int dir) const
    {
        f[0] = 0.5 * u[0] * u[0];
        return std::abs(u[0]);
    }
    //=============================================================================
    EulerFlux::EulerFlux(const EquationOfState *eos, int dim)
        : FluxFunction(dim, 2 + dim),
          eos_(eos)
    {
    }
    //=============================================================================
    const EquationOfState *EulerFlux::eos() const
    {
        return eos_;
    }
    //=============================================================================
    EulerFlux::State EulerFlux::get_state(std::span<const double> u) const
    {
        if ((int)u.size() != n_comp_)
        {
            // big problem !
        }

        State state{};
        state.rho = u[0];
        for (int i = 0; i < dim_; i++)
        {
            state.v[i] = u[i + 1] / state.rho;
            state.KE += state.v[i] * state.v[i];
        }
        state.KE *= 0.5 * state.rho;
        state.E = u[dim_ + 1];
        state.e = (state.E - state.KE) / state.rho;
        state.p = eos_->pressure(state.rho, state.e);

        return state;
    }
    //=============================================================================
    double EulerFlux::operator()(const State &state, std::span<double> f, int dir) const
    {
        f[0] = state.rho * state.v[dir];
        for (int i = 0; i < dim_; i++)
        {
            f[i + 1] = state.rho * state.v[i] * state.v[dir];
        }
        f[dir + 1] += state.p;
        f[dim_ + 1] = (state.E + state.p) * state.v[dir];

        return std::abs(state.v[dir]) + eos_->sound_speed(state.rho, state.p);
    }
    //=============================================================================
    double EulerFlux::operator()(std::span<const double> u, std::span<double> f, int dir) const
    {
        return operator()(get_state(u), f, dir);
    }
}