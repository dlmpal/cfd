#pragma once

#include "flux_function.hpp"
#include "grid.hpp"

namespace hc
{
    class NumericalFlux
    {
    public:
        NumericalFlux(const FluxFunction *flux);

        const FluxFunction *flux() const;

        virtual double compute_flux(std::span<const double> ul,
                                    std::span<const double> ur,
                                    std::span<double> f, int dir) const = 0;

    protected:
        const FluxFunction *flux_;
        mutable std::vector<double> fl_;
        mutable std::vector<double> fr_;
    };

    class LaxFriedrichsFlux : public NumericalFlux
    {
    public:
        LaxFriedrichsFlux(const FluxFunction *flux, const GridGeo *geo, const double *dt);

        double compute_flux(std::span<const double> ul,
                            std::span<const double> ur,
                            std::span<double> f, int dir) const override;

    private:
        const GridGeo *geo_;
        const double *dt_;
    };

    class ForceFlux : public NumericalFlux
    {
    public:
        ForceFlux(FluxFunction *flux, const GridGeo *geo, const double *dt);

        double compute_flux(std::span<const double> ul,
                            std::span<const double> ur,
                            std::span<double> f, int dir) const override;

    private:
        const GridGeo *geo_;
        const double *dt_;

        mutable std::vector<double> f_lf_;
        mutable std::vector<double> u_ri_;
        mutable std::vector<double> f_ri_;
    };

    class RusanovFlux : public NumericalFlux
    {
    public:
        RusanovFlux(FluxFunction *flux);

        double compute_flux(std::span<const double> ul,
                            std::span<const double> ur,
                            std::span<double> f, int dir) const override;
    };

    class HLLFlux : public NumericalFlux
    {
    public:
        HLLFlux(EulerFlux *flux);

        double compute_flux(std::span<const double> ul,
                            std::span<const double> ur,
                            std::span<double> f, int dir) const override;
    };

    class HLLCFlux : public NumericalFlux
    {
    public:
        HLLCFlux(EulerFlux *flux);

        double compute_flux(std::span<const double> ul,
                            std::span<const double> ur,
                            std::span<double> f, int dir) const override;

    private:
        mutable std::vector<double> u_hllc_l_;
        mutable std::vector<double> u_hllc_r_;
    };
}