#pragma once
#include <Eigen/Dense>
#include <memory>
#include "Constants.hpp"

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
		int BASIS_SIZE;
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

		// Computes the respective x or y component from a given input index
		inline int x(int idx) const {
			return idx / (2 * Constants::K_DISCRETIZATION) - Constants::K_DISCRETIZATION;
		};
		inline int y(int idx) const {
			return idx % (2 * Constants::K_DISCRETIZATION) - Constants::K_DISCRETIZATION;
		};
		inline int equal_up_to_Q(const Eigen::Vector2i& l, const Eigen::Vector2i& r) const {
			if (l == r) return 0;

			if (l(0) == r(0) + Constants::K_DISCRETIZATION || l(0) == r(0) - Constants::K_DISCRETIZATION) {
				if (l(1) == r(1) + Constants::K_DISCRETIZATION || l(1) == r(1) - Constants::K_DISCRETIZATION) {
					return 1;
				}
			}
			return -1;
		};
		inline void clean_factor_2pi(Eigen::Vector2i& toClean) const {
			// + Q is required for the modulo operation later
			// as well as referencing, which works on indizes from 0 to [2pi] and not from [-pi] to [pi]
			for (int i = 0; i < 2; i++)
			{
				toClean(i) += Constants::K_DISCRETIZATION;
				if(toClean(i) < 0){
					toClean(i) = ((2 * Constants::K_DISCRETIZATION) - abs(toClean(i) % (2 * Constants::K_DISCRETIZATION))) % (2 * Constants::K_DISCRETIZATION);
				} else {
					toClean(i) = (toClean(i) % (2 * Constants::K_DISCRETIZATION)) % (2 * Constants::K_DISCRETIZATION);
				}
			}
		};
		// Returns a c^+ c^+ (cc) type term, i.e. the SC or the eta order parameter
		inline double sc_type(Eigen::Vector2i left, Eigen::Vector2i right) const {
			right = -right;
			clean_factor_2pi(left);
			clean_factor_2pi(right);

			int offset = equal_up_to_Q(left, right);
			if (offset < 0) return 0;
			return expecs[2 + offset](right(0), right(1));
		};
		// Returns a c^+ c type term, i.e. the CDW order parameter or the number operator
		inline double cdw_type(Eigen::Vector2i left, Eigen::Vector2i right) const {
			clean_factor_2pi(left);
			clean_factor_2pi(right);

			int offset = equal_up_to_Q(left, right);
			if (offset < 0) return 0;		
			return expecs[offset](left(0), left(1));
		};


		inline double fermi_dirac(double energy) const {
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
		void computeCollectiveModes_v2(std::vector<std::vector<double>>& reciever);
	};
}