#define _USE_MATH_DEFINES

#include "UsingBroyden.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include "../Utility/Roots_Broyden.hpp"

using Eigen::MatrixXd;

namespace Hubbard {
	void UsingBroyden::data_set::print() const {
		std::cout << delta_cdw << "\t" << delta_sc << "\t" << delta_eta
			<< "\t" << sqrt(delta_cdw * delta_cdw + delta_sc * delta_sc + delta_eta * delta_eta) << std::endl;
	}

	double UsingBroyden::unperturbed_energy(double k_x, double k_y) const
	{
		return -2 * (cos(k_x) + cos(k_y));
	}

	void UsingBroyden::fillMatrix(double k_x, double k_y)
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

	UsingBroyden::UsingBroyden(ModelParameters& _params)
		: Model(_params), V(_params.V)
	{
		this->delta_cdw = abs(U - V) * 0.5;
		this->delta_sc = abs(U + V) * 0.5;
		if (V > 0) {
			this->delta_sc *= 0.25;
		}
		else {
			this->delta_cdw *= 0.25;
		}
		this->delta_eta = 0.01;

		this->hamilton = MatrixXd::Zero(4, 4);
	}

	UsingBroyden::data_set UsingBroyden::compute(const bool print/*=false*/)
	{
		MatrixXd rho = MatrixXd::Zero(4, 4);
		Eigen::SelfAdjointEigenSolver<MatrixXd> solver;

		auto lambda_func = [&](const Eigen::VectorXd& x, Eigen::VectorXd& F) {
			delta_cdw = x(0);
			delta_sc = x(1);
			delta_eta = x(2);
			F = Eigen::Vector3d::Zero();
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

					F(0) += rho(0, 1);
					F(1) += rho(0, 2);
					F(2) += rho(0, 3);
				}
			}
			setParameters(F);
			F -= x;
		};
		std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)> func = lambda_func;
		Eigen::VectorXd x0 = Eigen::Vector3d(delta_cdw, delta_sc, delta_eta);
		Eigen::VectorXd f0 = Eigen::Vector3d(delta_cdw, delta_sc, delta_eta);
		for (size_t i = 0; i < 30; i++)
		{
			func(x0, f0);
			x0(0) = delta_cdw;
			x0(1) = delta_sc;
			x0(2) = delta_eta;
		}

		if (!Utility::Roots::Broyden::compute(func, x0)) {
			std::cerr << "No convergence for [T, U, V] = [" << this->temperature << ", " << this->U << ", " << this->V << "]" << std::endl;
			delta_cdw = 0;
			delta_sc = 0;
			delta_eta = 0;
		}
		data_set ret = { this->delta_cdw, this->delta_sc, this->delta_eta };

		if (print) {
			func(x0, f0);
			std::cout << "x0 = (";
			for (int i = 0; i < 3; i++)
			{
				std::cout << " " << x0(i) << " ";
			}
			std::cout << ")\nf0 = (";
			for (int i = 0; i < 3; i++)
			{
				std::cout << " " << f0(i) << " ";
			}
			std::cout << ")   -> |f0| = " << f0.norm() << std::endl;
		}

		return ret;
	}
}