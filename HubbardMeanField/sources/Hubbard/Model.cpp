#define _USE_MATH_DEFINES

#include "Model.hpp"
#include "Constants.hpp"
#include <iostream>
#include <fstream>

namespace Hubbard {
	std::pair<int, int> Hubbard::Model::parseExpectationValue(std::string& str)
	{
		if (str[0] == 'n') {
			if (str[7] == '0') {
				return std::make_pair(0, 0);
			}
			else {
				return  std::make_pair(1, 1);
			}
		}
		if (str[0] == 'g') {
			if (str[3] == '0') {
				return std::make_pair(1, 0);
			}
			else {
				return std::make_pair(0, 1);
			}
		}
		if (str.substr(0, 3) == "eta") {
			str.erase(0, 4);
			if (str[0] == '1') {
				if (str[4] == '1') {
					return std::make_pair(0, 3);
				}
				else {
					return std::make_pair(1, 2);
				}
			}
			else {
				if (str[4] == '1') {
					return std::make_pair(3, 0);
				}
				else {
					return std::make_pair(2, 1);
				}
			}
		}
		if (str[0] == 'f') {
			str.erase(0, 2);
			if (str[0] == '1') {
				if (str[4] == '1') {
					return std::make_pair(1, 3);
				}
				else {
					return std::make_pair(0, 2);
				}
			}
			else {
				if (str[4] == '1') {
					return std::make_pair(3, 1);
				}
				else {
					return std::make_pair(2, 0);
				}
			}
		}
		return std::pair<int, int>(-1, -1);
	}

	void Model::parseCommutatorData()
	{
		struct Term {
			double prefactor = 0;
			std::pair<int, int> matrix_indizes;
			Term* next_term = nullptr;

			double computeValue(const Eigen::MatrixXd& expectation_values) const {
				double value = prefactor * expectation_values(matrix_indizes.first, matrix_indizes.second);
				if (next_term != nullptr) {
					value += next_term->computeValue(expectation_values);
				}
				return value;
			};
		};
		typedef std::pair<double, std::vector<Term>> coeff_term;
		std::vector<std::vector<coeff_term>> expressions;

		std::ifstream n_file("../data/commuting_N.txt");
		std::string line;
		while (std::getline(n_file, line)) {
			if (line[0] == ']') continue;
			else if (line[0] == '[') {
				// We get an entirely new expression
				expressions.push_back(std::vector<coeff_term>());
			}
			else if (line[0] == '}') continue;
			else if (line[0] == '{') {
				// We get a new coefficient within the old expression
				expressions.back().push_back(coeff_term());
			}
			else if (line[0] == '\t') {
				if (line[1] == '}') continue;
				else if (line[1] == '{') {
					// We get a new Wick-type term
					expressions.back().back().second.push_back(Term());
				}
				else if (line[1] == '\t') {
					if (line[2] == '{') {
						// Parse Wick
						line.erase(0, 3);
						if (line[0] == '+') {
						}
						else if (line[0] == '-') {
						}
						else {
							expressions.back().back().second.back().matrix_indizes = this->parseExpectationValue(line);
						}
					}
					else {
						// Parse prefactor
						expressions.back().back().second.back().prefactor = std::stod(line.erase(0, 2));
					}
				}
				else {
					// Parse coefficient
					int coeff = std::stoi(line.erase(0, 1));
					switch (coeff) {
					case -1:
						expressions.back().back().first = 1;
						break;
					case 0: // epsilon is k-dependend, hence we need further checks down the line
						expressions.back().back().first = -128;
						break;
					case 1:
						expressions.back().back().first = this->delta_cdw;
						break;
					case 2:
						expressions.back().back().first = this->delta_sc;
						break;
					case 3: // complex conjugate, but it's 0 for now anyways
						expressions.back().back().first = this->delta_eta;
						break;
					case 4:
						expressions.back().back().first = this->delta_eta;
						break;
					}
				}
			}
		}
	}

	Model::Model(double _temperature, double _U)
		: temperature(_temperature), U(_U)
	{
		this->delta_cdw = 0.1;
		this->delta_sc = 0.1;
		this->delta_eta = 0.001;
	}

	Model::Model(ModelParameters& _params)
		: temperature(_params.temperature), U(_params.U)
	{
		this->delta_cdw = 0.1;
		this->delta_sc = 0.1;
		this->delta_eta = 0.001;
	}

	void Model::getEnergies(std::vector<std::vector<double>>& reciever, double direction)
	{
		reciever.reserve(2 * Constants::K_DISCRETIZATION);
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
		double k_val = 0;
		for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
		{
			k_val = k * M_PI / Constants::K_DISCRETIZATION;
			fillMatrix(k_val, direction * k_val);
			solver.compute(hamilton, false);
			reciever.push_back(std::vector<double>(solver.eigenvalues().data(), solver.eigenvalues().data() + solver.eigenvalues().size()));
		}
	}

	Model::ModelParameters::ModelParameters(double _temperature, double _U, double _V, double _global_step, double _second_step,
		std::string _global_iterator_type, std::string _second_iterator_type)
		: global_iterator_type(_global_iterator_type), second_iterator_type(_second_iterator_type),
		global_step(_global_step), second_step(_second_step), temperature(_temperature), U(_U), V(_V)
	{
		if (second_iterator_type == "T") {
			second_it_min = temperature;
		}
		else if (second_iterator_type == "U") {
			second_it_min = U;
		}
		else if (second_iterator_type == "V") {
			second_it_min = V;
		}
		else {
			second_it_min = 0;
		}
	}

	// // // // // // // // // // // //
	void Model::ModelParameters::incrementer(std::string& s, const double step)
	{
		if (s == "T") {
			temperature += step;
		}
		else if (s == "U") {
			U += step;
		}
		else if (s == "V") {
			V += step;
		}
	}
	void Model::ModelParameters::setSecondIterator(int it_num)
	{
		if (second_iterator_type == "T") {
			temperature = second_it_min + it_num * second_step;
		}
		else if (second_iterator_type == "U") {
			U = second_it_min + it_num * second_step;
		}
		else if (second_iterator_type == "V") {
			V = second_it_min + it_num * second_step;
		}
	}
	void Model::ModelParameters::incrementGlobalIterator()
	{
		incrementer(global_iterator_type, global_step);
		setSecondIterator(0);
	}
	void Model::ModelParameters::incrementSecondIterator()
	{
		incrementer(second_iterator_type, second_step);
	}
	void Model::ModelParameters::printGlobal() const
	{
		std::cout << global_iterator_type << " = ";
		if (global_iterator_type == "T") {
			std::cout << temperature;
		}
		else if (global_iterator_type == "U") {
			std::cout << U;
		}
		else if (global_iterator_type == "V") {
			std::cout << V;
		}
	}
	void Model::data_set::print() const {
		std::cout << delta_cdw << "\t" << delta_sc << "\t" << delta_eta
			<< "\t" << sqrt(delta_cdw * delta_cdw + delta_sc * delta_sc + delta_eta * delta_eta) << std::endl;
	}
}