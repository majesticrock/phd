#pragma once
#include <Eigen/Dense>

namespace Hubbard {
	class Model
	{
	protected:
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
			void incrementGlobalIterator();
			void incrementSecondIterator();
			void printGlobal() const;
		};

		Model(double _temperature, double _U);
		Model(ModelParameters& _params);
	};
}