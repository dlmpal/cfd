#include "equation_of_state.hpp"
#include <cmath>

namespace hc
{
    //=============================================================================
    EquationOfState::EquationOfState()
    {
    }
    //=============================================================================
    IdealGasLaw::IdealGasLaw(double gamma)
        : gamma_(gamma)
    {
    }
    //=============================================================================
    double IdealGasLaw::gamma() const
    {
        return gamma_;
    }
    //=============================================================================
    double IdealGasLaw::sound_speed(double rho, double p) const
    {
        return std::sqrt(gamma_ * p / rho);
    }
    //=============================================================================
    double IdealGasLaw::pressure(double rho, double e) const
    {
        return (gamma_ - 1) * rho * e;
    }
}