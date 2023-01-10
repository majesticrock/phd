#define _USE_MATH_DEFINES

#include "UsingBroyden.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include "../Utility/Roots_Broyden.hpp"

using Eigen::MatrixXd;

namespace Hubbard {
	double UsingBroyden::unperturbed_energy(double k_x, double k_y) const
	{
		return -2 * (cos(k_x) + cos(k_y));
	}

	void UsingBroyden::fillHamiltonian(double k_x, double k_y)
	{
		hamilton.fill(0);

		Eigen::Vector4d cos_sin;
		cos_sin << cos(k_x), cos(k_y), sin(k_x), sin(k_y);
		hamilton(0, 1) = delta_cdw;
		hamilton(0, 2) = delta_sc + V * old_sc.dot(cos_sin);
		hamilton(0, 3) = delta_eta - V * old_eta.dot(cos_sin);

		hamilton(1, 2) = delta_eta + V * old_eta.dot(cos_sin);
		hamilton(1, 3) = delta_sc - V * old_sc.dot(cos_sin);
		hamilton(2, 3) = -delta_cdw;

		Eigen::MatrixXd buffer = hamilton.transpose();
		hamilton += buffer;
		double eps = unperturbed_energy(k_x, k_y);
		hamilton(0, 0) = eps;
		hamilton(1, 1) = -eps;
		hamilton(2, 2) = -eps;
		hamilton(3, 3) = eps;
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
		this->delta_eta = 0;

		this->hamilton = MatrixXd::Zero(4, 4);
		this->old_eta = Eigen::Vector4d::Zero();
		this->old_sc = Eigen::Vector4d::Zero();
	}

	Model::data_set UsingBroyden::computePhases(const bool print/*=false*/)
	{
		MatrixXd rho = MatrixXd::Zero(4, 4);
		Eigen::SelfAdjointEigenSolver<MatrixXd> solver;

		auto lambda_func = [&](const Eigen::VectorXd& x, Eigen::VectorXd& F) {
			delta_cdw = x(0);
			delta_sc = x(1);
			delta_eta = x(2);
			new_sc = Eigen::Vector4d::Zero();
			new_eta = Eigen::Vector4d::Zero();

			F = Eigen::Vector3d::Zero();
			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
				{
					double k_x = ((k + l) * (0.5 * M_PI)) / Constants::K_DISCRETIZATION;
					double k_y = ((k - l) * (0.5 * M_PI)) / Constants::K_DISCRETIZATION;
					fillHamiltonian(((k + l) * (0.5 * M_PI)) / Constants::K_DISCRETIZATION, ((k - l) * (0.5 * M_PI)) / Constants::K_DISCRETIZATION);
					solver.compute(hamilton);
					rho.fill(0);
					for (int i = 0; i < 4; i++)
					{
						rho(i, i) = fermi_dirac(solver.eigenvalues()(i));
					}
					rho = solver.eigenvectors() * rho * (solver.eigenvectors().transpose());
					F(0) += (rho(1, 0) - rho(2, 3));
					F(1) += (rho(0, 2) + rho(1, 3));
					F(2) += (rho(0, 3) + rho(1, 2));
					//std::cout << rho(0, 2) << "     " << rho(1, 3) << std::endl;
					new_sc(0) += cos(k_x) * (rho(0,2) - rho(1,3));
					new_sc(1) += cos(k_y) * (rho(0,2) - rho(1,3));
					new_sc(2) += sin(k_x) * (rho(0,2) - rho(1,3));
					new_sc(3) += sin(k_y) * (rho(0,2) - rho(1,3));

					new_eta(0) += cos(k_x) * (rho(0,3) - rho(1,2));
					new_eta(1) += cos(k_y) * (rho(0,3) - rho(1,2));
					new_eta(2) += sin(k_x) * (rho(0,3) - rho(1,2));
					new_eta(3) += sin(k_y) * (rho(0,3) - rho(1,2));
				}
			}
			old_sc = new_sc / (4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);
			old_eta = new_eta / (4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);
			setParameters(F);
			F -= x;
		};
		std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)> func = lambda_func;
		Eigen::VectorXd x0 = Eigen::Vector3d(delta_cdw, delta_sc, delta_eta);
		Eigen::VectorXd f0 = Eigen::Vector3d(delta_cdw, delta_sc, delta_eta);
		for (size_t i = 0; i < 200 && f0.squaredNorm() > 1e-15; i++)
		{
			func(x0, f0);
			if (abs(x0(0) + delta_cdw) < 1e-10) {
				delta_cdw = 0;
			}
			if (abs(x0(1) + delta_sc) < 1e-10) {
				delta_sc = 0;
			}
			if (abs(x0(2) + delta_eta) < 1e-10) {
				delta_eta = 0;
			}
			x0(0) = 0.5 * (delta_cdw + x0(0));
			x0(1) = 0.5 * (delta_sc + x0(1));
			x0(2) = 0.5 * (delta_eta + x0(2));
		}

		if (!Utility::Roots::Broyden::compute(func, x0)) {
			std::cerr << "No convergence for [T, U, V] = [" << this->temperature << ", " << this->U << ", " << this->V << "]" << std::endl;
			delta_cdw = 0;
			delta_sc = 0;
			delta_eta = 0;
		}
		data_set ret = { this->delta_cdw, this->delta_sc + old_sc.sum(), this->delta_eta + old_eta.sum()};
		std::cout << old_sc.sum() << std::endl;
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