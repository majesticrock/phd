#pragma once
#include <Eigen/Dense>
#include <memory>

namespace Hubbard {
	struct Term {
		double prefactor = 1;
		std::vector<std::pair<int, int>> matrix_indizes;
		std::vector<std::shared_ptr<Term>> next_terms;

		double computeValue(const Eigen::MatrixXd& expectation_values) const {
			double value = 0;
			double buf = 0;
			for (int i = 0; i < matrix_indizes.size(); i++)
			{
				if (matrix_indizes[i].first < 0) {
					value += 1;
				}
				else {
					buf = expectation_values(matrix_indizes[i].first, matrix_indizes[i].second);
					if (i < next_terms.size()) {
						buf *= next_terms[i]->computeValue(expectation_values);
					}
					value += buf;
				}
			}
			return prefactor * value;
		};
	};

	std::ostream& operator<<(std::ostream& os, const Term& t);

	class Model
	{
	private:
		std::pair<int, int> parseExpectationValue(std::string& str);
		void parseWick(std::shared_ptr<Term> lastOne, std::string& line);
	protected:
		int BASIS_SIZE = (2 * Constants::K_DISCRETIZATION) * (2 * Constants::K_DISCRETIZATION);
		typedef std::pair<double, std::vector<std::shared_ptr<Term>>> coeff_term;
		std::vector<std::vector<coeff_term>> expressions_M, expressions_N;

		double delta_sc, delta_cdw, delta_eta;
		Eigen::MatrixXd hamilton;
		double temperature;
		double U;

		/*
		* 0 - number operator
		* 1 - cdw
		* 2 - sc
		* 3 - eta
		*/
		std::vector<Eigen::MatrixXd> expecs;
		double sum_of_all[4] = { 0, 0, 0, 0 };
		// Quartic expecs - contracted using Wick's theorem (exact on MF level)
		std::vector<Eigen::MatrixXd> quartic;

		Eigen::MatrixXd M;
		Eigen::MatrixXd N;

		inline double fermi_dirac(double energy) {
			if (temperature > 1e-7) {
				return (1. / (1 + exp(energy / temperature)));
			}
			else {
				return ((energy > 0) ? 0 : 1);
			}
		};

		virtual double unperturbed_energy(double k_x, double k_y) const = 0;
		virtual void fillHamiltonian(double k_x, double k_y) = 0;

		virtual void compute_quartics();
		virtual void fill_M_N();
	public:
		class ModelParameters {
		private:
			std::string global_iterator_type;
			std::string second_iterator_type;
			double global_step;
			double second_step;
			double second_it_min;

			void incrementer(std::string& s, const double step);
		public:
			double temperature;
			double U;
			double V;

			ModelParameters(double _temperature, double _U, double _V, double global_step, double second_step,
				std::string _global_iterator_type, std::string _second_iterator_type);
			void setSecondIterator(int it_num);
			void incrementGlobalIterator();
			void incrementSecondIterator();
			inline double getGlobal() const {
				if (global_iterator_type == "T") {
					return temperature;
				}
				else if (global_iterator_type == "U") {
					return U;
				}
				else if (global_iterator_type == "V") {
					return V;
				}
				return -128;
			};
			void printGlobal() const;
		};
		struct data_set {
			double delta_cdw, delta_sc, delta_eta;
			void print() const;
		};

		Model(double _temperature, double _U);
		Model(ModelParameters& _params);
		// reciever is the vector the resulting data will be stored in
		// direction gives the angle between the k-path and the k_x axis in multiples of M_PI
		void getEnergies(std::vector<std::vector<double>>& reciever, double direction);
		virtual data_set computePhases(const bool print = false) = 0;
		void parseCommutatorData();
		void computeCollectiveModes(std::vector<std::vector<double>>& reciever, double direction);
		// version 2 use the non mean field hamilton for the commutation,
		// but the mean field system to obtain the expectation values
		void computeCollectiveModes_v2(std::vector<std::vector<double>>& reciever, double direction);
	};
}