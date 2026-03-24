#pragma once
#include "GlobalDefinitions.hpp"
#include "ModelAttributes.hpp"
#include <string>
#include <map>

#include <mrock/utility/InputFileReader.hpp>
#include <mrock/symbolic_operators/WickTerm.hpp>

namespace DWave {
    struct Model {
        typedef Eigen::VectorXd ParameterVector;

        const l_float dwave_coupling_in; ///< V_in
        const l_float phonon_coupling_in; ///< g_in
        const l_float fermi_energy; ///< in units of W
        const l_float omega_debye; ///< debye frequency
        const int N; ///< Number of discretization points
        
        const l_float dwave_coupling; ///< V_in / rho_F
        const l_float phonon_coupling; ///< g_in / rho_F
        l_float beta; ///< inverse temperature in units of W
        l_float chemical_potential; ///< chemical potential in units of W
        l_float filling_at_zero_temp; ///< filling at zero temperature

        ModelAttributes<l_float> Delta;

        Model(mrock::utility::InputFileReader& input);

        void iteration_step(const ParameterVector& initial_values, ParameterVector& result);

        inline int ravel_index(int kx, int ky) const {
            return kx + N * ky;
        }
        inline int unravel_x(int k) const {
            return k % N;
        }
        inline int unravel_y(int k) const {
            return k / N;
        }

        l_float dwave_factor(int kx, int ky) const;

        l_float dispersion(int kx, int ky) const;
        l_float quasiparticle_dispersion(int kx, int ky) const;

        l_float sc_expectation_value(int kx, int ky) const;
        l_float occupation_number(int kx, int ky) const;

        l_float compute_filling(const l_float mu) const;
        
        const std::map<mrock::symbolic_operators::OperatorType, std::vector<l_float>>& get_expectation_values() const;

        l_float delta_max() const;
        l_float delta_true() const;

        std::string info() const;
        std::string to_folder() const;

    private:
		mutable std::map<mrock::symbolic_operators::OperatorType, std::vector<l_float>> _expecs;
        const bool guaranteed_E_F{};
    };
}