#pragma once

#include "grid_array.hpp"
#include <functional>

namespace hc
{
    /// @brief Explicit Runge-Kutta integrator
    class ERKIntegrator
    {
    public:
        /// @brief Right-hand-side function F = F(t, S)
        using RHSFunction = std::function<void(const GridArray &S, GridArray &rhs, double time)>;

        /// @brief Create an Explicit Runge-Kutta integrator
        /// @param state State vector
        /// @param rhs RHS function
        /// @param n_stages Number of stages
        /// @param nodes Butcher tableau nodes
        /// @param weights Butcher tableau weights
        /// @param coeffs Butcher tableau coefficients
        ERKIntegrator(const GridArray &state, RHSFunction rhs,
                      int n_stages, std::vector<double> &&nodes,
                      std::vector<double> &&weights,
                      std::vector<std::vector<double>> &&coeffs);

        /// @brief Advance the state forward by one timestep
        /// @param S_old Old timestep state vector
        /// @param S_new New timestep state vector
        /// @param time Current time
        /// @param dt Current timestep size
        void advance(const GridArray &S_old, GridArray &S_new, double time, double dt);

    private:
        /// @brief RHS function
        RHSFunction rhs_;

        /// @brief Number of stages
        int n_stages_;

        /// @brief Butcher tableau nodes (ci)
        std::vector<double> nodes_;

        /// @brief Butcher tableau weights (bi)
        std::vector<double> weights_;

        /// @brief Butcher tableau coefficients (aij)
        std::vector<std::vector<double>> coeffs_;

        /// @brief RHS vector for the intermediate stages
        std::vector<GridArray> stages_;
    };

    /// @brief Available Explicit Runge-Kutta integrator types
    enum class ERKType
    {
        fe,  // Forward Euler
        rk2, // Second-order Runge-Kutta
        rk3, // Third-order Runge-Kutta
        rk4  // Classic fourth-order Runge-Kutta
    };

    /// @brief Create an ERK integrator of specific type
    ERKIntegrator create_erk(const GridArray &state, ERKIntegrator::RHSFunction rhs, ERKType type);
}