#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <iomanip>
#include "BasicHubbardModel.hpp"

using Eigen::MatrixXd;

namespace Hubbard {
	double BasicHubbardModel::unperturbed_energy(double k_x, double k_y) const
	{
		return -2 * (cos(k_x) + cos(k_y));
	}

	void BasicHubbardModel::fillMatrix(double k_x, double k_y)
	{
		hamilton.fill(0);

		hamilton(0, 1) = delta_cdw;
		hamilton(0, 2) = delta_sc;
		hamilton(0, 3) = delta_eta;
		hamilton(1, 2) = delta_eta;
		hamilton(1, 3) = delta_sc;
		hamilton(2, 3) = -delta_cdw;

		Eigen::MatrixXd buffer = hamilton.transpose();
		hamilton += buffer;

		hamilton(0, 0) = unperturbed_energy(k_x, k_y);
		hamilton(1, 1) = unperturbed_energy(k_x + M_PI, k_y + M_PI);
		hamilton(2, 2) = -unperturbed_energy(-k_x, -k_y);
		hamilton(3, 3) = -unperturbed_energy(-k_x + M_PI, -k_y + M_PI);
	}

	BasicHubbardModel::BasicHubbardModel(ModelParameters& _params)
		: Model(_params)
	{
		this->delta_sc = 0.1;
		this->delta_cdw = 0.1;
		this->delta_eta = 0.01;

		this->hamilton = MatrixXd::Zero(4, 4);
	}

	Hubbard::Model::data_set BasicHubbardModel::compute(const bool print/*=false*/)
	{
		MatrixXd rho = MatrixXd::Zero(4, 4);
		Eigen::SelfAdjointEigenSolver<MatrixXd> solver;
		const double EPSILON = 1e-8;
		double sc = 0, cdw = 0, eta = 0;
		double old_parameters[3] = { 100, 100, 100 };
		double error = 100;
		double error_cdw = 100, error_cdw_osc = 100;
		double error_sc = 100, error_sc_osc = 100;
		double error_eta = 100, error_eta_osc = 100;
		const int MAX_STEPS = 10000;

		for (size_t i = 0; i < MAX_STEPS && error > EPSILON; i++)
		{
			sc = 0;
			cdw = 0;
			eta = 0;

			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
				{
					fillMatrix((k * M_PI) / Constants::K_DISCRETIZATION, (l * M_PI) / Constants::K_DISCRETIZATION);
					solver.compute(hamilton);
					rho.fill(0);

					for (int i = 0; i < 4; i++)
					{
						rho(i, i) = fermi_dirac(solver.eigenvalues()(i));
					}
					rho = solver.eigenvectors() * rho * (solver.eigenvectors().transpose());

					cdw += rho(0, 1);
					sc += rho(0, 2);
					eta += rho(0, 3);
				}
			}

			old_parameters[0] = delta_cdw;
			old_parameters[1] = delta_sc;
			old_parameters[2] = delta_eta;
			setParameters(cdw, sc, eta);
			error_cdw = abs(delta_cdw - old_parameters[0]);
			error_sc = abs(delta_sc - old_parameters[1]);
			error_eta = abs(delta_eta - old_parameters[2]);
			error = error_cdw + error_sc + error_eta;

			delta_cdw = 0.5 * old_parameters[0] + 0.5 * delta_cdw;
			delta_sc = 0.5 * old_parameters[1] + 0.5 * delta_sc;
			delta_eta = 0.5 * old_parameters[2] + 0.5 * delta_eta;

			error_cdw_osc = abs(delta_cdw + old_parameters[0]);
			error_sc_osc = abs(delta_sc + old_parameters[1]);
			error_eta_osc = abs(delta_eta + old_parameters[2]);

			delta_cdw = ((error_cdw_osc) < EPSILON) ? 0 : delta_cdw;
			delta_sc = ((error_sc_osc) < EPSILON) ? 0 : delta_sc;
			delta_eta = ((error_eta_osc) < EPSILON) ? 0 : delta_eta;

			if (print) {
				double total = 0;
				for (size_t i = 0; i < 3; i++)
				{
					total += old_parameters[i] * old_parameters[i];
				}
				std::cout << i << ":\t" << std::fixed << std::setprecision(8)
					<< delta_cdw << "\t" << delta_sc << "\t" << delta_eta << "\t" << sqrt(total) << "\t" << error << std::endl;
			}
			if (i == MAX_STEPS - 1) {
				std::cerr << "[T, U] = [" << this->temperature << ", " << this->U << "]\tConvergence at " << error << std::endl;
				delta_cdw = 0;
				delta_sc = 0;
				delta_eta = 0;
			}
		}
		data_set ret;
		ret.delta_cdw = delta_cdw;
		ret.delta_sc = delta_sc;
		ret.delta_eta = delta_eta;

		return ret;
	}
}