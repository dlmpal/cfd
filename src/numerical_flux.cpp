#include "numerical_flux.hpp"
#include "equation_of_state.hpp"

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
          f_lf_(flux_->n_comp()),
          u_ri_(flux_->n_comp()),
          f_ri_(flux_->n_comp())
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
            f_lf_[i] = 0.5 * (fr_[i] + fl_[i]) - 0.5 * geo_->dX(dir) / dt * (ur[i] - ul[i]);
        }

        // Richtmeyer flux
        for (int i = 0; i < flux_->n_comp(); i++)
        {
            u_ri_[i] = 0.5 * (ur[i] + ul[i]) - 0.5 * dt / geo_->dX(dir) * (fr_[i] - fl_[i]);
        }
        const double v_ri = (*flux_)(u_ri_, f_ri_, dir);

        for (int i = 0; i < flux_->n_comp(); i++)
        {
            f[i] = 0.5 * (f_lf_[i] + f_ri_[i]);
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
    //=============================================================================
    HLLFlux::HLLFlux(EulerFlux *flux)
        : NumericalFlux(flux)
    {
    }
    //=============================================================================
    double HLLFlux::compute_flux(std::span<const double> ul,
                                 std::span<const double> ur,
                                 std::span<double> f, int dir) const
    {
        const EulerFlux *euler = static_cast<const EulerFlux *>(flux_);
        const EquationOfState *eos = euler->eos();
        const int dim = flux_->dim();

        const EulerFlux::State wl = euler->get_state(ul);
        const EulerFlux::State wr = euler->get_state(ur);

        (*euler)(wl, fl_, dir);
        (*euler)(wr, fr_, dir);

        const double cl = eos->sound_speed(wl.rho, wl.p);
        const double cr = eos->sound_speed(wr.rho, wr.p);

        const double s_plus = std::max(std::abs(wl.v[dir]) + cl, std::abs(wr.v[dir]) + cr);
        const double sl = -s_plus;
        const double sr = s_plus;

        for (int i = 0; i < flux_->n_comp(); i++)
        {
            if (sl > 0)
            {
                f[i] = fl_[i];
            }
            else if (sr > 0)
            {
                f[i] = (sr * fl_[i] - sl * fr_[i] + sl * sr * (ur[i] - ul[i])) / (sr - sl);
            }
            else
            {
                f[i] = fr_[i];
            }
        }

        return std::max(std::abs(sl), std::abs(sr));
    }
    //=============================================================================
    HLLCFlux::HLLCFlux(EulerFlux *flux)
        : NumericalFlux(flux),
          u_hllc_l_(flux->n_comp()),
          u_hllc_r_(flux->n_comp())
    {
    }
    //=============================================================================
    double HLLCFlux::compute_flux(std::span<const double> ul,
                                  std::span<const double> ur,
                                  std::span<double> f, int dir) const
    {
        const EulerFlux *euler = static_cast<const EulerFlux *>(flux_);
        const EquationOfState *eos = euler->eos();
        const int dim = flux_->dim();

        const EulerFlux::State wl = euler->get_state(ul);
        const EulerFlux::State wr = euler->get_state(ur);

        (*euler)(wl, fl_, dir);
        (*euler)(wr, fr_, dir);

        const double cl = eos->sound_speed(wl.rho, wl.p);
        const double cr = eos->sound_speed(wr.rho, wr.p);

        const double s_plus = std::max(std::abs(wl.v[dir]) + cl, std::abs(wr.v[dir]) + cr);
        const double sl = -s_plus;
        const double sr = s_plus;

        const double s = (wr.p - wl.p + wl.rho * wl.v[dir] * (sl - wl.v[dir]) - wr.rho * wr.v[dir] * (sr - wr.v[dir])) /
                         (wl.rho * (sl - wr.v[dir]) - wr.rho * (sr - wr.v[dir]));

        u_hllc_l_[0] = wl.rho * (sl - wl.v[dir]) / (sl - s);
        u_hllc_r_[0] = wr.rho * (sr - wr.v[dir]) / (sr - s);
        u_hllc_r_[dim + 1] = u_hllc_r_[0] * (wr.E / wr.rho + (s - wr.v[dir]) * (s + wr.p / wr.rho / (sr - wr.v[dir])));
        u_hllc_l_[dim + 1] = u_hllc_l_[0] * (wl.E / wl.rho + (s - wl.v[dir]) * (s + wl.p / wl.rho / (sl - wl.v[dir])));
        for (int i = 0; i < dim; i++)
        {
            if (i == dir)
            {
                u_hllc_l_[i + 1] = u_hllc_l_[0] * s;
                u_hllc_r_[i + 1] = u_hllc_r_[0] * s;
            }
            else
            {
                u_hllc_l_[i + 1] = u_hllc_l_[0] * wl.v[i];
                u_hllc_r_[i + 1] = u_hllc_r_[0] * wr.v[i];
            }
        }

        for (int i = 0; i < flux_->n_comp(); i++)
        {
            if (sl > 0)
            {
                f[i] = fl_[i];
            }
            else if (s > 0)
            {
                f[i] = fl_[i] + sl * (u_hllc_l_[i] - ul[i]);
            }
            else if (sr > 0)
            {
                f[i] = fr_[i] + sr * (u_hllc_r_[i] - ur[i]);
            }
            else
            {
                f[i] = fr_[i];
            }
        }

        return std::max(std::abs(sl), std::abs(sr));
    }
}