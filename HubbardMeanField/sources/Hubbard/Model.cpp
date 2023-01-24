#define _USE_MATH_DEFINES

#include "Model.hpp"
#include "Constants.hpp"
#include <iostream>
#include <fstream>

using Eigen::MatrixXd;

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

	Model::Model(double _temperature, double _U)
		: temperature(_temperature), U(_U)
	{
		this->delta_cdw = 0.1;
		this->delta_sc = 0.1;
		this->delta_eta = 0.001;
	}

	Model::Model(ModelParameters& _params)
		: temperature(_params.temperature), U(_params.U)
	{
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

	void Model::computeCollectiveModes(std::vector<std::vector<double>>& reciever, double direction)
	{
		parseCommutatorData();
		reciever.reserve(2 * Constants::K_DISCRETIZATION);
		const int BASIS_SIZE = 16;

		Eigen::MatrixXd M = Eigen::MatrixXd::Zero(BASIS_SIZE, BASIS_SIZE);
		Eigen::MatrixXd N = Eigen::MatrixXd::Zero(BASIS_SIZE, BASIS_SIZE);
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

				if (++j >= BASIS_SIZE) j = ++i;
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

	void Model::computeCollectiveModes_v2(std::vector<std::vector<double>>& reciever, double direction)
	{
		// First off we need to compute every possible expectation value
		// We use the mean field system's symmetries
		// i.e. there are only the standard SC, CDW, Eta and N operators non-zero
		// spin symmetry and conservation of momentum up to addition of Q holds
		MatrixXd rho = MatrixXd::Zero(4, 4);
		Eigen::SelfAdjointEigenSolver<MatrixXd> solver;
		/*
		* 0 - number operator
		* 1 - cdw
		* 2 - sc
		* 3 - eta
		*/
		std::vector<MatrixXd> expecs(4, MatrixXd::Zero(2 * Constants::K_DISCRETIZATION, 2 * Constants::K_DISCRETIZATION));
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
				}
			}
		}
		int BASIS_SIZE = (2 * Constants::K_DISCRETIZATION) * (2 * Constants::K_DISCRETIZATION);
		// Quartic expecs - contracted using Wick's theorem (exact on MF level)
		std::vector<MatrixXd> quartic(14, MatrixXd::Zero(BASIS_SIZE, BASIS_SIZE));

		// Computes the respective x or y component from a given input index
		auto x = [&](int idx) -> int {
			return idx / (2 * Constants::K_DISCRETIZATION) - Constants::K_DISCRETIZATION;
		};
		auto y = [&](int idx) -> int {
			return idx % (2 * Constants::K_DISCRETIZATION) - Constants::K_DISCRETIZATION;
		};
		auto equal_up_to_Q = [&](const std::pair<int, int>& l, const std::pair<int, int>& r) -> int {
			if (l == r) return 0;

			if (l.first == r.first + Constants::K_DISCRETIZATION || l.first == r.first - Constants::K_DISCRETIZATION) {
				if (l.second == r.second + Constants::K_DISCRETIZATION || l.second == r.second - Constants::K_DISCRETIZATION) {
					return 1;
				}
			}
			return -1;
		};

		// Returns a c^+ c^+ (cc) type term, i.e. the SC or the eta order parameter
		auto sc_type = [&](const std::pair<int, int>& left, const std::pair<int, int>& right) -> double {
			int offset = equal_up_to_Q(left, right);
			if (offset < 0) return 0;
			return expecs[2 + offset](right.first, right.second);
		};
		// Returns a c^+ c type term, i.e. the CDW order parameter or the number operator
		auto cdw_type = [&](const std::pair<int, int>& left, const std::pair<int, int>& right) -> double {
			int offset = equal_up_to_Q(left, right);
			if (offset < 0) return 0;
			return expecs[offset](right.first, right.second);
		};

		for (int k = 0; k < BASIS_SIZE; k++)
		{
			std::pair<int, int> pair_k = { x(k), y(k) };
			for (int l = 0; l < BASIS_SIZE; l++)
			{
				std::pair<int, int> pair_l = { x(l), y(l) };
				std::pair<int, int> pair_p = { 0, 0 };

				std::pair<int, int> left, right, left2, right2;
				quartic[0](k, l) += sc_type(left, right) * sc_type(left2, right2);
			}
		}
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