#include "Model.hpp"

#include <mrock/utility/better_to_string.hpp>
#include <boost/math/tools/roots.hpp>

#include <limits>

namespace DWave {
    Model::Model(mrock::utility::InputFileReader& input) 
        : 
        dwave_coupling_in{ input.getDouble("dwave_coupling") },
        phonon_coupling_in{ input.getDouble("phonon_coupling") },
        fermi_energy{ input.getDouble("fermi_energy") },
        omega_debye{ input.getDouble("omega_debye") },
        N{input.getInt("N")},
        dwave_coupling{  dwave_coupling_in  / (N*N) },
        phonon_coupling{ phonon_coupling_in / (N*N) },
        beta{ input.getDouble("beta") },
        chemical_potential{ fermi_energy },
        filling_at_zero_temp{ 0. },
        Delta{ decltype(Delta)::Random(N*N + 1) }
    {
        Delta[N*N] = fermi_energy;

        for (int kx = 0; kx < N; ++kx) {
            for (int ky = 0; ky < N; ++ky) {
                //std::cout << kx << "\t" << ky << "\t" << dispersion(kx, ky) << "\t" << fermi_function(dispersion(kx, ky), -1) << std::endl;
                filling_at_zero_temp += fermi_function(dispersion(kx, ky), -1);
            }
        }
        filling_at_zero_temp /= (N*N);
        std::cout << "Filling at zero temperature: " << filling_at_zero_temp << std::endl;
        std::cout << "Compare initial sc-filling:  " << compute_filling(fermi_energy) << std::endl;
    }

    void Model::iteration_step(const ParameterVector& initial_values, ParameterVector& result)
    {
        static int step_num = 0;
        result.setZero();
        this->Delta.fill_with(initial_values);
        this->chemical_potential = initial_values(N*N);
        this->get_expectation_values();

        auto fit_occupation = [&](const l_float mu) -> l_float {
            return filling_at_zero_temp - compute_filling(mu);
        };
        auto bracket_root = [&]() -> std::pair<l_float, l_float> {
            constexpr l_float step = 0.005;
            l_float a{chemical_potential - step};
            l_float b{chemical_potential + step};

            l_float fa{fit_occupation(a)};
            l_float fb{fit_occupation(b)};

            l_float dfa, dfb;

            if (fa * fb < 0.0) {
                return std::pair<l_float, l_float>{a, b};
            }
            else {
                dfa = fit_occupation(a - step) - fa;
                dfb = fit_occupation(b + step) - fb;
            }

            while (fa * fb > 0.0) {
                if (fb * dfb < 0.0){
                    b += step;
                    fb = fit_occupation(b);
                    dfb = fit_occupation(b + step) - fb;
                }
                if (fa * dfa < 0.0) {
                    a -= step;
                    fa = fit_occupation(a);
                    dfa = fit_occupation(a - step) - fa;
                }
            }
            return std::pair<l_float, l_float>{a, b};
        };

        std::uintmax_t boost_max_it{100U};
        const auto bracket = bracket_root();
        const auto best_mu = boost::math::tools::toms748_solve(fit_occupation, 
                bracket.first, bracket.second,
                boost::math::tools::eps_tolerance<l_float>(32), boost_max_it);

        const l_float fa{fit_occupation(best_mu.first)};
        const l_float fb{fit_occupation(best_mu.second)};
        result(N*N) = is_zero(fb - fa) ? 0.5 * (best_mu.first + best_mu.second) : (best_mu.first * fb - best_mu.second * fa) / (fb - fa);
        

#pragma omp parallel for
        for (int k = 0; k < N*N; k++)
        {
            const l_float energy_k = dispersion(unravel_x(k), unravel_y(k));
            l_float __part{};
            for (int l = 0; l < N*N; ++l)
            {
                if (std::abs(energy_k - dispersion(unravel_x(l), unravel_y(l))) <= omega_debye)
                {
                    __part -= _expecs[mrock::symbolic_operators::SC_Type][l];
                }
            }
            result(k) = phonon_coupling * __part;
            __part = l_float{};
            for (int l = 0; l < N*N; ++l)
            {
                __part -= _expecs[mrock::symbolic_operators::SC_Type][l] * dwave_factor(unravel_x(l), unravel_y(l));
            }
            result(k) += dwave_coupling * dwave_factor(unravel_x(k), unravel_y(k)) * __part;
        }
        this->Delta.fill_with(result, 0.5);
        
        this->chemical_potential = this->Delta[N*N];
		this->Delta.clear_noise(PRECISION);
		result -= initial_values;
		++step_num;
    }

    l_float Model::dwave_factor(int kx, int ky) const
    {
        return (std::cos((pi * (2 * kx - 1)) / N) - std::cos((pi * (2 * ky - 1)) / N));
    }

    l_float Model::dispersion(int kx, int ky) const
    {
        return -0.5 * (std::cos((pi * (2 * kx - 1)) / N) + std::cos((pi * (2 * ky - 1)) / N)) - chemical_potential;
    }

    l_float Model::quasiparticle_dispersion(int kx, int ky) const
    {
        const l_float eps = dispersion(kx, ky);
        return std::sqrt(eps*eps + Delta[ravel_index(kx, ky)] * Delta[ravel_index(kx, ky)]);
    }


    l_float Model::sc_expectation_value(int kx, int ky) const
    {
        if (is_zero(Delta[ravel_index(kx, ky)])) return l_float{};
        if (beta < 0)
		    return -Delta[ravel_index(kx, ky)] / (2 * quasiparticle_dispersion(kx, ky));
        else
            return -Delta[ravel_index(kx, ky)] / (2 * quasiparticle_dispersion(kx, ky)) * std::tanh(0.5 * beta * quasiparticle_dispersion(kx, ky));
    }

    l_float Model::occupation_number(int kx, int ky) const
    {
        const l_float eps_mu = dispersion(kx, ky);
		if (is_zero(Delta[ravel_index(kx, ky)])) {
			return fermi_function(eps_mu, beta);
		}
        const l_float E = quasiparticle_dispersion(kx, ky);
        if (beta < 0)
            return 0.5 * (1 - (eps_mu / E));
        else
		    return 0.5 * (1 - (eps_mu / E) * std::tanh(0.5 * beta * E));
    }

    l_float Model::compute_filling(const l_float mu) const
    {
        auto number_operator = [&](int kx, int ky) -> l_float {
            const l_float eps = -0.5 * (std::cos((pi * (2 * kx - 1)) / N) + std::cos((pi * (2 * ky - 1)) / N)) - mu;
            const l_float E = sqrt(eps*eps + Delta[ravel_index(kx, ky)]*Delta[ravel_index(kx, ky)]);
            if (is_zero(Delta[ravel_index(kx, ky)]))
			    return fermi_function(eps, beta);
            if (beta < 0)
                return 0.5 * (1 - (eps / E));
            else
		        return 0.5 * (1 - (eps / E) * std::tanh(0.5 * beta * E));
        };
        l_float n{};
        for (int kx = 0; kx < N; ++kx) {
            for (int ky = 0; ky < N; ++ky) {
                n += number_operator(kx, ky);
            }
        }
        n /= (N*N);
        return n;
    }


    const std::map<mrock::symbolic_operators::OperatorType, std::vector<l_float>> &Model::get_expectation_values() const
    {
		if (_expecs.empty()) {
			_expecs.emplace(mrock::symbolic_operators::Number_Type, std::vector<l_float>(N*N));
			_expecs.emplace(mrock::symbolic_operators::SC_Type, std::vector<l_float>(N*N));
		}
#pragma omp parallel for
		for (int k = 0; k < N*N; ++k) {
			_expecs.at(mrock::symbolic_operators::Number_Type)[k] = this->occupation_number(unravel_x(k), unravel_y(k));
			_expecs.at(mrock::symbolic_operators::SC_Type)[k]  = this->sc_expectation_value(unravel_x(k), unravel_y(k));
		}

		return _expecs;
	}

    l_float Model::delta_max() const
    {
        return std::abs(*std::max_element(Delta.begin(), Delta.begin() + N*N,
			[](const l_float& lhs, const l_float& rhs) {
				return std::abs(lhs) < std::abs(rhs);
			}
		));
    }

    l_float Model::delta_true() const
    {
        l_float current_min = 1000.;
        for (int kx = 0; kx < N; ++kx) {
            for (int ky = 0; ky < N; ++ky) {
                const l_float E = quasiparticle_dispersion(kx, ky);
                if (E < current_min) {
                    current_min = E;
                }
            }
        }
        return current_min;
    }

    std::string Model::info() const
    {
        return "Model: g=" + std::to_string(phonon_coupling) + "\tV=" + std::to_string(dwave_coupling_in);
    }

    std::string Model::to_folder() const
    {
        auto improved_string = [](l_float number) -> std::string {
			if (std::floor(number) == number) {
				// If the number is a whole number, format it with one decimal place
				std::ostringstream out;
				out.precision(1);
				out << std::fixed << number;
				return out.str();
			}
			else {
				std::string str = mrock::utility::better_to_string(number, std::chars_format::fixed);
				// Remove trailing zeroes
				str.erase(str.find_last_not_of('0') + 1, std::string::npos);
				str.erase(str.find_last_not_of('.') + 1, std::string::npos);
				return str;
			}
			};
        
        return "/N=" + std::to_string(N) + "/"
            + "V=" + improved_string(dwave_coupling_in) + "/"
            + "E_F=" + improved_string(fermi_energy) + "/"
            + "g=" + improved_string(phonon_coupling_in) + "/"
            + "omega_D=" + improved_string(omega_debye) + "/";
    }
}