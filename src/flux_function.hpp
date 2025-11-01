#pragma once

#include <array>
#include <cmath>
#include <span>
#include <vector>

namespace hc
{
    class FluxFunction
    {
    public:
        FluxFunction(int dim, int n_comp);

        int dim() const;

        int n_comp() const;

        virtual double operator()(std::span<const double> u, std::span<double> f, int dir) const = 0;

    protected:
        int dim_;
        int n_comp_;
    };

    /// @todo Extend to 2D and 3D
    class AdvectionFlux : public FluxFunction
    {
    public:
        AdvectionFlux(double speed);

        double speed() const;

        double operator()(std::span<const double> u, std::span<double> f, int dir) const override;

    private:
        double speed_;
    };

    /// @todo Extend to 2D and 3D
    class BurgersFlux : public FluxFunction
    {
    public:
        BurgersFlux();

        double operator()(std::span<const double> u, std::span<double> f, int dir) const override;
    };

    // Forward declartion
    class EquationOfState;

    class EulerFlux : public FluxFunction
    {
    public:
        struct State
        {
            double rho;
            std::array<double, 3> v;
            double KE;
            double E;
            double e;
            double p;
        };

        EulerFlux(const EquationOfState *eos, int dim);

        const EquationOfState *eos() const;

        State get_state(std::span<const double> U) const;

        double operator()(const State &state, std::span<double> f, int dir) const;

        double operator()(std::span<const double> u, std::span<double> f, int dir) const override;

    private:
        const EquationOfState *eos_;
    };
}