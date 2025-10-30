#pragma once

namespace hc
{
    class EquationOfState
    {
    public:
        EquationOfState();

        virtual double gamma() const = 0;
        virtual double sound_speed(double rho, double p) const = 0;
        virtual double pressure(double rho, double e) const = 0;
    };

    class IdealGasLaw : public EquationOfState
    {
    public:
        IdealGasLaw(double gamma);

        double gamma() const override;
        double sound_speed(double rho, double p) const override;
        double pressure(double rho, double e) const override;

    private:
        /// @brief Adiabatic index
        double gamma_;
    };
}