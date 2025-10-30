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
    double EulerFlux::operator()(std::span<const double> u, std::span<double> f, int dir) const
    {
        if ((int)u.size() != n_comp_)
        {
            // big problem !
        }

        // Get state variables
        const double rho = u[0];
        std::array<double, 3> v;
        double KE = 0.0;
        for (int i = 0; i < dim_; i++)
        {
            v[i] = u[i + 1] / rho;
            KE += 0.5 * rho * v[i] * v[i];
        }
        const double E = u[dim_ + 1];
        const double e = (E - KE) / rho;
        const double p = eos_->pressure(rho, e);

        // Compute fluxes
        f[0] = rho * v[dir];
        for (int i = 0; i < dim_; i++)
        {
            f[i + 1] = rho * v[i] * v[dir];
        }
        f[dir + 1] += p;
        f[dim_ + 1] = (E + p) * v[dir];

        const double speed = std::sqrt(2.0 / rho * KE); /// @todo Check validity
        const double sound = eos_->sound_speed(rho, p);

        return speed + sound;
    }
}