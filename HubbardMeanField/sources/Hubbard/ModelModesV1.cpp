#include "Model.hpp"
#include <iostream>
#include <fstream>

namespace Hubbard {
	std::pair<int, int> Hubbard::Model::parseExpectationValue(std::string& str)
	{
		if (str[0] == 'n') {
			if (str[7] == '0') {
				return std::make_pair(0, 0);
			}
			else {
				return  std::make_pair(1, 1);
			}
		}
		if (str[0] == 'g') {
			if (str[3] == '0') {
				return std::make_pair(1, 0);
			}
			else {
				return std::make_pair(0, 1);
			}
		}
		if (str.substr(0, 3) == "eta") {
			str.erase(0, 4);
			if (str[0] == '1') {
				if (str[4] == '1') {
					return std::make_pair(0, 3);
				}
				else {
					return std::make_pair(1, 2);
				}
			}
			else {
				if (str[4] == '1') {
					return std::make_pair(3, 0);
				}
				else {
					return std::make_pair(2, 1);
				}
			}
		}
		if (str[0] == 'f') {
			str.erase(0, 2);
			if (str[0] == '1') {
				if (str[4] == '1') {
					return std::make_pair(1, 3);
				}
				else {
					return std::make_pair(0, 2);
				}
			}
			else {
				if (str[4] == '1') {
					return std::make_pair(3, 1);
				}
				else {
					return std::make_pair(2, 0);
				}
			}
		}
		return std::pair<int, int>(-1, -1);
	}

	void Hubbard::Model::parseWick(std::shared_ptr<Term> lastOne, std::string& line)
	{
		if (line[0] == '+' || line[0] == '-') {
			int pre = (line[0] == '+') ? 1 : -1;
			line.erase(0, 3);
			lastOne->matrix_indizes.push_back(this->parseExpectationValue(line));
			//line.erase(0, line.find('{') + 1);

			lastOne->next_terms.push_back(std::make_shared<Term>(Term()));
			lastOne->next_terms.back()->prefactor = pre;
			line.erase(0, line.find('{') + 1);
			parseWick(lastOne->next_terms.back(), line);
			line.erase(0, line.find('}') + 1);

			while (line[0] == ',') {
				pre = (line[2] == '+') ? 1 : -1;
				line.erase(0, 5);
				lastOne->matrix_indizes.push_back(this->parseExpectationValue(line));
				lastOne->next_terms.push_back(std::make_shared<Term>(Term()));
				lastOne->next_terms.back()->prefactor = pre;

				line.erase(0, line.find('{') + 1);
				parseWick(lastOne->next_terms.back(), line);
				line.erase(0, line.find('}') + 1);
			}
		}
		else {
			lastOne->matrix_indizes.push_back(this->parseExpectationValue(line));
			line.erase(0, line.find('}') + 1);
		}
	}

	void Model::parseCommutatorData()
	{
		std::shared_ptr<Term> last_term;

		std::string line;
		std::ifstream m_file("../data/commuting_M.txt");
		while (std::getline(m_file, line)) {
			if (line[0] == ']') continue;
			else if (line[0] == '[') {
				// We get an entirely new expression
				expressions_M.push_back(std::vector<coeff_term>());
			}
			else if (line[0] == '}') continue;
			else if (line[0] == '{') {
				// We get a new coefficient within the old expression
				expressions_M.back().push_back(coeff_term());
			}
			else if (line[0] == '\t') {
				if (line[1] == '}') continue;
				else if (line[1] == '{') {
					// We get a new Wick-type term
					expressions_M.back().back().second.push_back(std::make_shared<Term>(Term()));
					last_term = expressions_M.back().back().second.back();
				}
				else if (line[1] == '\t') {
					if (line[2] == '{') {
						// Parse Wick
						line.erase(0, 3);
						this->parseWick(last_term, line);
					}
					else if (line[2] == 'I') {
						last_term->matrix_indizes.push_back(std::make_pair<int, int>(-1, -1));
					}
					else {
						// Parse prefactor
						last_term->prefactor = std::stod(line.erase(0, 2));
					}
				}
				else {
					// Parse coefficient
					int coeff = std::stoi(line.erase(0, 1));
					switch (coeff) {
					case -1:
						expressions_M.back().back().first = 1;
						break;
					case 0: // epsilon is k-dependend, hence we need further checks down the line
						expressions_M.back().back().first = -128;
						break;
					case 1:
						expressions_M.back().back().first = this->delta_cdw;
						break;
					case 2:
						expressions_M.back().back().first = this->delta_sc;
						break;
					case 3: // complex conjugate, but it's 0 for now anyways
						expressions_M.back().back().first = this->delta_eta;
						break;
					case 4:
						expressions_M.back().back().first = this->delta_eta;
						break;
					}
				}
			}
		}

		std::ifstream n_file("../data/commuting_N.txt");
		while (std::getline(n_file, line)) {
			if (line[0] == ']') continue;
			else if (line[0] == '[') {
				// We get an entirely new expression
				expressions_N.push_back(std::vector<coeff_term>());
			}
			else if (line[0] == '}') continue;
			else if (line[0] == '{') {
				// We get a new coefficient within the old expression
				expressions_N.back().push_back(coeff_term());
			}
			else if (line[0] == '\t') {
				if (line[1] == '}') continue;
				else if (line[1] == '{') {
					// We get a new Wick-type term
					expressions_N.back().back().second.push_back(std::make_shared<Term>(Term()));
					last_term = expressions_N.back().back().second.back();
				}
				else if (line[1] == '\t') {
					if (line[2] == '{') {
						// Parse Wick
						line.erase(0, 3);
						this->parseWick(last_term, line);
					}
					else if (line[2] == 'I') {
						last_term->matrix_indizes.push_back(std::make_pair<int, int>(-1, -1));
					}
					else {
						// Parse prefactor
						last_term->prefactor = std::stod(line.erase(0, 2));
					}
				}
				else {
					// Parse coefficient
					int coeff = std::stoi(line.erase(0, 1));
					switch (coeff) {
					case -1:
						expressions_N.back().back().first = 1;
						break;
					case 0: // epsilon is k-dependend, hence we need further checks down the line
						expressions_N.back().back().first = -128;
						break;
					case 1:
						expressions_N.back().back().first = this->delta_cdw;
						break;
					case 2:
						expressions_N.back().back().first = this->delta_sc;
						break;
					case 3: // complex conjugate, but it's 0 for now anyways
						expressions_N.back().back().first = this->delta_eta;
						break;
					case 4:
						expressions_N.back().back().first = this->delta_eta;
						break;
					}
				}
			}
		}
	}

	void Model::computeCollectiveModes(std::vector<std::vector<double>>& reciever, double direction)
	{
		parseCommutatorData();
		reciever.reserve(2 * Constants::K_DISCRETIZATION);
		const int BASIS_SIZE_LOCAL = 16;

		Eigen::MatrixXd M = Eigen::MatrixXd::Zero(BASIS_SIZE_LOCAL, BASIS_SIZE_LOCAL);
		Eigen::MatrixXd N = Eigen::MatrixXd::Zero(BASIS_SIZE_LOCAL, BASIS_SIZE_LOCAL);
		Eigen::Matrix4d rho = Eigen::Matrix4d::Zero();

		auto fill_matrices = [&](double kx, double ky) {
			int i = 0;
			int j = 0;
			for (int n = 0; n < expressions_M.size(); n++)
			{
				double value = 0;
				for (const auto& b : expressions_M[n]) {
					double buffer = 0;
					for (const auto& c : b.second) {
						buffer += c->computeValue(rho);
					}
					if (b.first != -128) {
						value += b.first * buffer;
					}
					else {
						value += unperturbed_energy(kx, ky) * buffer;
					}
				}
				if (abs(value) < 1e-12) {
					value = 0;
				}
				M(i, j) = value;
				M(j, i) = value;

				value = 0;
				for (const auto& b : expressions_N[n]) {
					double buffer = 0;
					for (const auto& c : b.second) {
						buffer += c->computeValue(rho);
					}
					if (b.first != -128) {
						value += b.first * buffer;
					}
					else {
						value += unperturbed_energy(kx, ky) * buffer;
					}
				}
				if (abs(value) < 1e-12) {
					value = 0;
				}
				N(i, j) = value;
				N(j, i) = value;

				if (++j >= BASIS_SIZE_LOCAL) j = ++i;
			}
		};

		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
		Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> gen_solver;
		Eigen::EigenSolver<Eigen::MatrixXd> solver_2;
		for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
		{
			double k_x = cos(M_PI * direction) * k * M_PI / Constants::K_DISCRETIZATION;
			double k_y = sin(M_PI * direction) * k * M_PI / Constants::K_DISCRETIZATION;
			fillHamiltonian(k_x, k_y);
			solver.compute(hamilton);
			rho.fill(0);
			for (int i = 0; i < 4; i++)
			{
				rho(i, i) = fermi_dirac(solver.eigenvalues()(i));
			}
			rho = solver.eigenvectors() * rho * (solver.eigenvectors().transpose());
			fill_matrices(k_x, k_y);

			//Eigen::MatrixXd b = M;
			//M = N;
			//N = b;

			Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> COD(N);
			Eigen::MatrixXd c = COD.pseudoInverse() * M;

			N += 1e-10 * Eigen::MatrixXd::Identity(N.rows(), N.cols());
			//M += 1e-10 * Eigen::MatrixXd::Identity(N.rows(), N.cols());
			solver.compute(N, false);
			for (int i = 0; i < solver.eigenvalues().size(); i++)
			{
				if (solver.eigenvalues()(i) < 0) {
					//std::cerr << "k=(" << k_x / M_PI << "," << k_y / M_PI << ")   N is not positive!\n" << N << std::endl << std::endl;
					break;
				}
			}
			gen_solver.compute(M, N);
			solver_2.compute(c);

			if (k == -Constants::K_DISCRETIZATION + 100) {
				std::cout << hamilton << std::endl << std::endl;
				std::cout << "##############################\n\n";
				std::cout << M << "\n\n" << N << "\n\n" << c << std::endl << std::endl;
				std::cout << gen_solver.eigenvalues() << "\n\n" << gen_solver.eigenvectors() << std::endl;
			}
			//reciever.push_back(std::vector<double>(gen_solver.eigenvalues().data(), gen_solver.eigenvalues().data() + gen_solver.eigenvalues().size()));
			Eigen::VectorXd ev = solver_2.eigenvalues().real();
			if (solver_2.eigenvalues().imag().squaredNorm() > 1e-12) {
				std::cout << solver_2.eigenvalues() << std::endl;
			}
			reciever.push_back(std::vector<double>(ev.data(), ev.data() + ev.size()));

			//std::sort(reciever.back().begin(), reciever.back().end());
		}
	}
}