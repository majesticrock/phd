#pragma once
#include <Eigen/Dense>

namespace Hubbard {
	class Model
	{
	protected:
		double delta_sc, delta_cdw, delta_eta;
		Eigen::MatrixXd hamilton;
		double temperature;
		double U;

		inline double fermi_dirac(double energy) {
			if (temperature > 1e-7) {
				return (1. / (1 + exp(energy / temperature)));
			}
			else {
				return ((energy > 0) ? 0 : 1);
			}
		};

		virtual double unperturbed_energy(double k_x, double k_y) const = 0;
		virtual void fillMatrix(double k_x, double k_y) = 0;
		std::pair<int, int> parseExpectationValue(std::string& str);
		void parseCommutatorData();
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
			void printGlobal() const;
		};
		struct data_set {
			double delta_cdw, delta_sc, delta_eta;
			void print() const;
		};

		Model(double _temperature, double _U);
		Model(ModelParameters& _params);
		// reciever is the vector the resulting data will be stored in
		// direction gives how strongly k_y is varied compared to x
		// i.e. 0 means k_y does not vary at all while 1 means k_x and k_y vary at the same rate
		// which results in a 45° angle in the first BZ
		void getEnergies(std::vector<std::vector<double>>& reciever, double direction);
	};
}