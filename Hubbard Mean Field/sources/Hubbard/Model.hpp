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
		Model(double _temperature, double _U);
	};
}