#include "numerical_flux.hpp"
#include <iostream>

namespace hc
{
    //=============================================================================
    NumericalFlux::NumericalFlux(const FluxFunction *flux)
        : flux_(flux),
          fl_(flux_->n_comp()),
          fr_(flux_->n_comp())
    {
    }
    //=============================================================================
    const FluxFunction *NumericalFlux::flux() const
    {
        return flux_;
    }
    //=============================================================================
    LaxFriedrichsFlux::LaxFriedrichsFlux(const FluxFunction *flux, const GridGeo *geo, const double *dt)
        : NumericalFlux(flux),
          geo_(geo),
          dt_(dt)
    {
    }
    //=============================================================================
    double LaxFriedrichsFlux::compute_flux(std::span<const double> ul,
                                           std::span<const double> ur,
                                           std::span<double> f, int dir) const
    {
        const auto sl = (*flux_)(ul, fl_, dir);
        const auto sr = (*flux_)(ur, fr_, dir);

        for (int i = 0; i < flux_->n_comp(); i++)
        {
            f[i] = 0.5 * (fr_[i] + fl_[i]) - 0.5 * geo_->dX(dir) / *dt_ * (ur[i] - ul[i]);
        }

        return std::max(sl, sr);
    }
    //=============================================================================
    ForceFlux::ForceFlux(FluxFunction *flux, const GridGeo *geo, const double *dt)
        : NumericalFlux(flux),
          geo_(geo),
          dt_(dt),
          f_lf(flux_->n_comp()),
          u_ri(flux_->n_comp()),
          f_ri(flux_->n_comp())
    {
    }
    //=============================================================================
    double ForceFlux::compute_flux(std::span<const double> ul,
                                   std::span<const double> ur,
                                   std::span<double> f, int dir) const
    {
        const double sl = (*flux_)(ul, fl_, dir);
        const double sr = (*flux_)(ur, fr_, dir);
        const double dt = *dt_;

        // Lax-Friedrichs flux
        for (int i = 0; i < flux_->n_comp(); i++)
        {
            f_lf[i] = 0.5 * (fr_[i] + fl_[i]) - 0.5 * geo_->dX(dir) / dt * (ur[i] - ul[i]);
        }

        // Richtmeyer flux
        for (int i = 0; i < flux_->n_comp(); i++)
        {
            u_ri[i] = 0.5 * (ur[i] + ul[i]) - 0.5 * dt / geo_->dX(dir) * (fr_[i] - fl_[i]);
        }
        const double v_ri = (*flux_)(u_ri, f_ri, dir);

        for (int i = 0; i < flux_->n_comp(); i++)
        {
            f[i] = 0.5 * (f_lf[i] + f_ri[i]);
        }

        return std::max(sl, sr);
    }
    //=============================================================================
    RusanovFlux::RusanovFlux(FluxFunction *flux)
        : NumericalFlux(flux)
    {
    }
    //=============================================================================
    double RusanovFlux::compute_flux(std::span<const double> ul,
                                     std::span<const double> ur,
                                     std::span<double> f, int dir) const
    {
        const double sl = (*flux_)(ul, fl_, dir);
        const double sr = (*flux_)(ur, fr_, dir);

        const double smax = std::max(sl, sr);

        for (int i = 0; i < flux_->n_comp(); i++)
        {
            f[i] = 0.5 * (fr_[i] + fl_[i]) - 0.5 * smax * (ur[i] - ul[i]);
        }

        return smax;
    }
}