#define _USE_MATH_DEFINES
#define _NUM(momentum) (expecs[0](x(momentum) + Constants::K_DISCRETIZATION, y(momentum) + Constants::K_DISCRETIZATION))
#define _CDW(momentum) (expecs[1](x(momentum) + Constants::K_DISCRETIZATION, y(momentum) + Constants::K_DISCRETIZATION))
#define _SC(momentum) (expecs[2](x(momentum) + Constants::K_DISCRETIZATION, y(momentum) + Constants::K_DISCRETIZATION))
#define _ETA(momentum) (expecs[3](x(momentum) + Constants::K_DISCRETIZATION, y(momentum) + Constants::K_DISCRETIZATION))

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
		hamilton(0, 2) = delta_sc;//+ V * old_sc.dot(cos_sin);
		hamilton(0, 3) = delta_eta;// - V * old_eta.dot(cos_sin);

		hamilton(1, 2) = delta_eta;// + V * old_eta.dot(cos_sin);
		hamilton(1, 3) = delta_sc;// -V * old_sc.dot(cos_sin);
		hamilton(2, 3) = -delta_cdw;

		Eigen::MatrixXd buffer = hamilton.transpose();
		hamilton += buffer;
		double eps = unperturbed_energy(k_x, k_y);
		hamilton(0, 0) = eps;
		hamilton(1, 1) = -eps;
		hamilton(2, 2) = -eps;
		hamilton(3, 3) = eps;
	}

	void UsingBroyden::compute_quartics()
	{
		Model::compute_quartics();

		Eigen::Vector2i vec_k, vec_l, vec_Q;
		vec_Q = { Constants::K_DISCRETIZATION, Constants::K_DISCRETIZATION };

		// Computes the respective x or y component from a given input index
		auto x = [&](int idx) -> int {
			return idx / (2 * Constants::K_DISCRETIZATION) - Constants::K_DISCRETIZATION;
		};
		auto y = [&](int idx) -> int {
			return idx % (2 * Constants::K_DISCRETIZATION) - Constants::K_DISCRETIZATION;
		};
		auto equal_up_to_Q = [&](const Eigen::Vector2i& l, const Eigen::Vector2i& r) -> int {
			if (l == r) return 0;

			if (l(0) == r(0) + Constants::K_DISCRETIZATION || l(0) == r(0) - Constants::K_DISCRETIZATION) {
				if (l(1) == r(1) + Constants::K_DISCRETIZATION || l(1) == r(1) - Constants::K_DISCRETIZATION) {
					return 1;
				}
			}
			return -1;
		};// Returns a c^+ c^+ (cc) type term, i.e. the SC or the eta order parameter
		auto sc_type = [&](Eigen::Vector2i left, Eigen::Vector2i right) -> double {
			int offset = equal_up_to_Q(left, right);
			if (offset < 0) return 0;
			right += vec_Q;
			while (right(0) < 0 || right(1) < 0) {
				right += 2 * vec_Q;
			}
			right(0) %= (2 * Constants::K_DISCRETIZATION);
			right(1) %= (2 * Constants::K_DISCRETIZATION);
			return expecs[2 + offset](right(0), right(1));
		};
		// Returns a c^+ c type term, i.e. the CDW order parameter or the number operator
		auto cdw_type = [&](Eigen::Vector2i left, Eigen::Vector2i right) -> double {
			int offset = equal_up_to_Q(left, right);
			if (offset < 0) return 0;
			left += vec_Q;
			while (left(0) < 0 || left(1) < 0) {
				left += 2 * vec_Q;
			}
			left(0) %= (2 * Constants::K_DISCRETIZATION);
			left(1) %= (2 * Constants::K_DISCRETIZATION);
			return expecs[offset](left(0), left(1));
		};

		auto f = [&](int kx, int ky) -> double {
			return (V / BASIS_SIZE) * 2 * (cos((M_PI * kx) / Constants::K_DISCRETIZATION) + cos((M_PI * ky) / Constants::K_DISCRETIZATION));
		};

		quartic.push_back(Eigen::MatrixXd::Zero(BASIS_SIZE, BASIS_SIZE));
		quartic.push_back(Eigen::MatrixXd::Zero(BASIS_SIZE, BASIS_SIZE));

		for (int k = 0; k < BASIS_SIZE; k++)
		{
			double buffer[4] = { 0,0,0,0 };
			for (int l = 0; l < BASIS_SIZE; l++)
			{
				vec_l = { x(l), y(l) };
				buffer[0] += f(vec_l(0), vec_l(1)) * sc_type(-vec_k - vec_l, vec_k + vec_l);
				buffer[1] += f(vec_l(0), vec_l(1)) * sc_type(-vec_k - vec_l, vec_k + vec_l + vec_Q);
				buffer[2] += f(vec_l(0), vec_l(1)) * cdw_type(-vec_k - vec_l, -vec_k - vec_l);
				buffer[3] += f(vec_l(0), vec_l(1)) * cdw_type(-vec_k - vec_l, -vec_k - vec_l + vec_Q);

				// number 6 -> 1a
				quartic[6](k, l) += sc_type(vec_l, -vec_l) * sc_type(-vec_k, vec_k);
				quartic[6](k, l) -= sc_type(vec_l, -vec_l - vec_Q) * sc_type(-vec_k, vec_k - vec_Q);
				quartic[6](k, l) += f(x(k) - x(l), y(k) - y(l)) * cdw_type(vec_l, vec_l) * cdw_type(-vec_k, -vec_k);
				quartic[6](k, l) += f(x(k) - x(l) + vec_Q(0), y(k) - y(l) + vec_Q(1)) * cdw_type(vec_l, vec_l - vec_Q) * cdw_type(-vec_k - vec_Q, -vec_k);
			}
			quartic[5](k, k) += sc_type(vec_k, -vec_k) * buffer[0];
			quartic[5](k, k) += sc_type(vec_k + vec_Q, -vec_k) * buffer[1];
			quartic[5](k, k) -= cdw_type(-vec_k, -vec_k) * buffer[2];
			quartic[5](k, k) -= cdw_type(-vec_k + vec_Q, -vec_k) * buffer[3];
			quartic[5](k, k) += 2 * cdw_type(-vec_k, -vec_k) * sum_of_all[0]; // cos(0) = 1
			quartic[5](k, k) -= 2 * cdw_type(-vec_k + vec_Q, -vec_k) * sum_of_all[1]; // minus because cos(pi) = -1
		}
	}

	void UsingBroyden::fill_M_N()
	{
		Model::fill_M_N();
		int BASIS_SIZE = (2 * Constants::K_DISCRETIZATION) * (2 * Constants::K_DISCRETIZATION);
		auto x = [&](int idx) -> int {
			return idx / (2 * Constants::K_DISCRETIZATION) - Constants::K_DISCRETIZATION;
		};
		auto y = [&](int idx) -> int {
			return idx % (2 * Constants::K_DISCRETIZATION) - Constants::K_DISCRETIZATION;
		};
		auto f = [&](int kx, int ky) -> double {
			return (V / BASIS_SIZE) * 2 * (cos((M_PI * kx) / Constants::K_DISCRETIZATION) + cos((M_PI * ky) / Constants::K_DISCRETIZATION));
		};

		for (int k = 0; k < BASIS_SIZE; k++)
		{
			M(k, k) += 1;
			for (int l = 0; l < BASIS_SIZE; l++)
			{
				M(l, k) += 2 * f(x(l) - x(k), y(l) - y(k)) * (1 - 2 * (_NUM(l) + _NUM(k) - quartic[0](l, k)));
			}
		}
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
		this->old_eta << 0.1, 0.1, 0.1, 0.1;
		this->old_sc << 0.1, 0.1, 0.1, 0.1;
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
					fillHamiltonian(k_x, k_y);
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

					//new_sc(0) += cos(k_x) * (rho(0, 2) - rho(1, 3));
					//new_sc(1) += cos(k_y) * (rho(0, 2) - rho(1, 3));
					//new_sc(2) += sin(k_x) * (rho(0, 2) - rho(1, 3));
					//new_sc(3) += sin(k_y) * (rho(0, 2) - rho(1, 3));
					//
					//new_eta(0) += cos(k_x) * (rho(0, 3) - rho(1, 2));
					//new_eta(1) += cos(k_y) * (rho(0, 3) - rho(1, 2));
					//new_eta(2) += sin(k_x) * (rho(0, 3) - rho(1, 2));
					//new_eta(3) += sin(k_y) * (rho(0, 3) - rho(1, 2));
				}
			}
			//old_sc = new_sc / (4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);
			//old_eta = new_eta / (4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);
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