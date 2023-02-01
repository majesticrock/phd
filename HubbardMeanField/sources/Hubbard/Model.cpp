#define _USE_MATH_DEFINES
#define _NUM(momentum) (expecs[0](x(momentum) + Constants::K_DISCRETIZATION, y(momentum) + Constants::K_DISCRETIZATION))

#include "Model.hpp"
#include <iostream>
#include <chrono>

namespace Hubbard {
	void Model::compute_quartics()
	{
		Eigen::Vector2i vec_k, vec_l;
		const Eigen::Vector2i vec_Q = { Constants::K_DISCRETIZATION, Constants::K_DISCRETIZATION };

		for (int k = 0; k < BASIS_SIZE; k++)
		{
			vec_k = { x(k), y(k) };
			for (int l = 0; l < BASIS_SIZE; l++)
			{
				vec_l = { x(l), y(l) };
				quartic[0](k, l) += sc_type(vec_l, -vec_l) * sc_type(-vec_k, vec_k);
				quartic[0](k, l) += sc_type(vec_l + vec_Q, -vec_l) * sc_type(-vec_k, vec_k + vec_Q);
				quartic[0](k, l) -= cdw_type(-vec_k, -vec_k) * cdw_type(-vec_l, -vec_l);
				quartic[0](k, l) -= cdw_type(-vec_k + vec_Q, -vec_k) * cdw_type(-vec_l, -vec_l + vec_Q);
				if (k == l) { // Factor 2 because of the spin sum
					quartic[0](k, l) += 2 * sum_of_all[0] * cdw_type(-vec_l, -vec_k);
				}
				else if (vec_k == vec_l + vec_Q || vec_k == vec_l - vec_Q) {
					quartic[0](k, l) += 2 * sum_of_all[1] * cdw_type(-vec_l, -vec_k);;
				}

				quartic[1](k, l) += sc_type(vec_l, -vec_l) * sc_type(-vec_k, vec_k);
				quartic[1](k, l) += sc_type(vec_l, -vec_l - vec_Q) * sc_type(-vec_k, vec_k - vec_Q);
				quartic[1](k, l) += cdw_type(vec_l, vec_l) * cdw_type(-vec_k, -vec_k);
				quartic[1](k, l) += cdw_type(vec_l, vec_l - vec_Q) * cdw_type(-vec_k - vec_Q, -vec_k);

				

				quartic[3](k, l) += cdw_type(-vec_k, -vec_k) * cdw_type(-vec_l, -vec_l);
				quartic[3](k, l) += cdw_type(-vec_k - vec_Q, -vec_k) * cdw_type(-vec_l, -vec_l - vec_Q);
				if (k == l) {
					quartic[3](k, l) -= sum_of_all[0] * cdw_type(vec_l, vec_k);
				}
				else if (vec_k == vec_l + vec_Q || vec_k == vec_l - vec_Q) {
					quartic[3](k, l) -= sum_of_all[1] * cdw_type(vec_l, vec_k);
				}
			}
			// the formulas include delta_kl, hence we dont need to recompute the same term for each l
			quartic[2](k, k) += sc_type(vec_k, -vec_k) * sum_of_all[2];
			quartic[2](k, k) += sc_type(vec_k, -vec_k + vec_Q) * sum_of_all[3];
			quartic[2](k, k) += cdw_type(vec_k, vec_k) * sum_of_all[0];
			quartic[2](k, k) += cdw_type(vec_k + vec_Q, vec_k) * sum_of_all[1];
		}

		//for(int c = 0; c < 4; c++) {
		//	Eigen::MatrixXd test = quartic[c] - quartic[c].transpose();
		//	for (size_t i = 0; i < test.rows(); i++)
		//	{
		//		vec_k = {x(i), y(i)};
		//		for (size_t j = i; j < test.cols(); j++)
		//		{
		//			vec_l = {x(j), y(j)};
		//			if (abs(test(i, j)) > 1e-8) {
		//				std::cout << c << "\t" << test(i, j) << "\t\t" << i << "\t" << j << std::endl;
		//				
		//				break;
		//			}
		//		}
		//	}
		//}
	}

	void Model::fill_M_N()
	{
		M = Eigen::MatrixXd::Zero(BASIS_SIZE, BASIS_SIZE);
		N = Eigen::MatrixXd::Zero(BASIS_SIZE, BASIS_SIZE);

		for (int k = 0; k < BASIS_SIZE; k++)
		{
			N(k, k) = -(1 - 2 * _NUM(k));

			M(k, k) += 2 * unperturbed_energy((M_PI * x(k)) / Constants::K_DISCRETIZATION, (M_PI * y(k)) / Constants::K_DISCRETIZATION)
				* (1 - 2 * _NUM(k));

			M(k, k) += (U / BASIS_SIZE) * (2 * sum_of_all[0] - quartic[2](k, k));

			for (int l = 0; l < BASIS_SIZE; l++)
			{
				M(l, k) += (U / BASIS_SIZE) * (1 - 2 * (_NUM(l) + _NUM(k) - quartic[1](l, k) - quartic[3](l, k)));
			}
		}
	}

	Model::Model(double _temperature, double _U)
		: temperature(_temperature), U(_U)
	{
		this->BASIS_SIZE = (2 * Constants::K_DISCRETIZATION) * (2 * Constants::K_DISCRETIZATION);
		this->delta_cdw = 0.1;
		this->delta_sc = 0.1;
		this->delta_eta = 0.001;
	}

	Model::Model(ModelParameters& _params)
		: temperature(_params.temperature), U(_params.U)
	{
		this->BASIS_SIZE = (2 * Constants::K_DISCRETIZATION) * (2 * Constants::K_DISCRETIZATION);
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
			fillHamiltonian(cos(M_PI * direction) * k_val, sin(M_PI * direction) * k_val);
			solver.compute(hamilton, false);
			reciever.push_back(std::vector<double>(solver.eigenvalues().data(), solver.eigenvalues().data() + solver.eigenvalues().size()));
		}
	}

	void Model::computeCollectiveModes_v2(std::vector<std::vector<double>>& reciever)
	{
		// First off we need to compute every possible expectation value
		// We use the mean field system's symmetries
		// i.e. there are only the standard SC, CDW, Eta and N operators non-zero
		// spin symmetry and conservation of momentum up to addition of Q holds
		std::chrono::time_point begin = std::chrono::steady_clock::now();
		std::chrono::time_point end = std::chrono::steady_clock::now();

		Eigen::MatrixXd rho = Eigen::MatrixXd::Zero(4, 4);
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;

		expecs = std::vector<Eigen::MatrixXd>(4, Eigen::MatrixXd::Zero(2 * Constants::K_DISCRETIZATION, 2 * Constants::K_DISCRETIZATION));
		quartic = std::vector<Eigen::MatrixXd>(4, Eigen::MatrixXd::Zero(BASIS_SIZE, BASIS_SIZE));

		for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
		{
			for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
			{
				fillHamiltonian((k * M_PI) / Constants::K_DISCRETIZATION, (l * M_PI) / Constants::K_DISCRETIZATION);
				solver.compute(hamilton);
				rho.fill(0);

				for (int i = 0; i < 4; i++)
				{
					rho(i, i) = fermi_dirac(solver.eigenvalues()(i));
				}
				rho = solver.eigenvectors() * rho * (solver.eigenvectors().transpose());
				for (int idx = 0; idx < 4; idx++)
				{
					expecs[idx](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = rho(idx, 0);
					sum_of_all[idx] += rho(idx, 0);
				}
			}
		}

		end = std::chrono::steady_clock::now();
		std::cout << "Time for expectation values: " 
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
		begin = std::chrono::steady_clock::now();

		compute_quartics();

		end = std::chrono::steady_clock::now();
		std::cout << "Time for quartic values: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
		begin = std::chrono::steady_clock::now();

		fill_M_N();

		end = std::chrono::steady_clock::now();
		std::cout << "Time for filling of M and N: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
		begin = std::chrono::steady_clock::now();

		const Eigen::MatrixXd test = M;
		for (size_t i = 0; i < test.rows(); i++)
		{
			for (size_t j = i + 1; j < test.cols(); j++)
			{
				if (abs(test(i, j)) > 1e-8) {
					std::cout << test(i, j) << "\t\t" << i << "\t" << j << std::endl;
				}
			}

			
		}

		Eigen::EigenSolver<Eigen::MatrixXd> gen_solver;
		Eigen::MatrixXd toSolve = M.inverse() * N;

		gen_solver.compute(toSolve, false);
		reciever.resize(1);
		Eigen::VectorXd ev = gen_solver.eigenvalues().real();
		reciever[0] = std::vector<double>(ev.data(), ev.data() + ev.size());

		for (size_t i = 0; i < ev.size(); i++)
		{
			if (abs(gen_solver.eigenvalues()(i).imag()) > 1e-6) {
				std::cout << gen_solver.eigenvalues()(i) << std::endl;
			}
		}

		std::sort(reciever.back().begin(), reciever.back().end());

		end = std::chrono::steady_clock::now();
		std::cout << "Time for solving M and N: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
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
		std::cout << global_iterator_type << " = " << getGlobal();
	}
	void Model::data_set::print() const {
		std::cout << delta_cdw << "\t" << delta_sc << "\t" << delta_eta
			<< "\t" << sqrt(delta_cdw * delta_cdw + delta_sc * delta_sc + delta_eta * delta_eta) << std::endl;
	}

	std::ostream& operator<<(std::ostream& os, const Term& t)
	{
		os << t.prefactor << " * [ ";
		for (int i = 0; i < t.matrix_indizes.size(); i++)
		{
			os << "(" << t.matrix_indizes[i].first << "," << t.matrix_indizes[i].second << ") ";
			if (i < t.next_terms.size()) {
				os << "* " << *(t.next_terms[i]);
			}
		}
		os << " ] ";
		return os;
	}
}