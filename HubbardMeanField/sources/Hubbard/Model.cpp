#define _USE_MATH_DEFINES
#define _NUM(momentum) (expecs[0](x(momentum) + Constants::K_DISCRETIZATION, y(momentum) + Constants::K_DISCRETIZATION))

#include "Model.hpp"
#include "Constants.hpp"
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

	void Model::compute_quartics()
	{
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
		};

		// Returns a c^+ c^+ (cc) type term, i.e. the SC or the eta order parameter
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

				quartic[2](k, l) += sc_type(vec_l, -vec_l) * sc_type(-vec_k, vec_k);
				quartic[2](k, l) += sc_type(vec_l, -vec_l - vec_Q) * sc_type(-vec_k - vec_Q, vec_k);
				if (k == l) {
					quartic[2](k, l) += sum_of_all[0] * cdw_type(vec_k, vec_k);
				}
				else if (vec_k == vec_l + vec_Q || vec_k == vec_l - vec_Q) {
					quartic[2](k, l) += sum_of_all[1] * cdw_type(vec_k, vec_k);
				}

				quartic[3](k, l) += cdw_type(-vec_k, -vec_k) * cdw_type(-vec_l, -vec_l);
				quartic[3](k, l) += cdw_type(-vec_k - vec_Q, -vec_k) * cdw_type(-vec_l, -vec_l - vec_Q);
				if (k == l) {
					quartic[3](k, l) -= sum_of_all[0] * cdw_type(vec_k, vec_k);
				}
				else if (vec_k == vec_l + vec_Q || vec_k == vec_l - vec_Q) {
					quartic[3](k, l) -= sum_of_all[1] * cdw_type(vec_k, vec_k);
				}
			}
			// the formulas include delta_kl, hence we dont need to recompute the same term for each l
			quartic[4](k, k) += sc_type(vec_k, -vec_k) * sum_of_all[2];
			quartic[4](k, k) += sc_type(vec_k, -vec_k + vec_Q) * sum_of_all[3];
			quartic[4](k, k) += cdw_type(vec_k, vec_k) * sum_of_all[0];
			quartic[4](k, k) += cdw_type(vec_k + vec_Q, vec_k) * sum_of_all[1];
		}
	}

	void Model::fill_M_N()
	{
		auto x = [&](int idx) -> int {
			return idx / (2 * Constants::K_DISCRETIZATION) - Constants::K_DISCRETIZATION;
		};
		auto y = [&](int idx) -> int {
			return idx % (2 * Constants::K_DISCRETIZATION) - Constants::K_DISCRETIZATION;
		};
		M = Eigen::MatrixXd::Zero(BASIS_SIZE, BASIS_SIZE);
		N = Eigen::MatrixXd::Zero(BASIS_SIZE, BASIS_SIZE);

		for (int k = 0; k < BASIS_SIZE; k++)
		{
			N(k, k) = -(1 - 2 * _NUM(k));

			M(k, k) = 2 * unperturbed_energy((M_PI * x(k)) / Constants::K_DISCRETIZATION,
				(M_PI * y(k)) / Constants::K_DISCRETIZATION) * (1 - 2 * _NUM(k));
			M(k, k) += (U / BASIS_SIZE) * (2 * sum_of_all[0] - quartic[4](k, k));

			for (int l = 0; l < BASIS_SIZE; l++)
			{
				M(l, k) += (U / BASIS_SIZE) * (1 - 2 * (_NUM(l) + _NUM(k) - quartic[1](l, k) - quartic[3](l, k)));
			}
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

	void Model::computeCollectiveModes_v2(std::vector<std::vector<double>>& reciever, double direction)
	{
		// First off we need to compute every possible expectation value
		// We use the mean field system's symmetries
		// i.e. there are only the standard SC, CDW, Eta and N operators non-zero
		// spin symmetry and conservation of momentum up to addition of Q holds
		Eigen::MatrixXd rho = Eigen::MatrixXd::Zero(4, 4);
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;

		expecs = std::vector<Eigen::MatrixXd>(4, Eigen::MatrixXd::Zero(2 * Constants::K_DISCRETIZATION, 2 * Constants::K_DISCRETIZATION));
		quartic = std::vector<Eigen::MatrixXd>(6, Eigen::MatrixXd::Zero(BASIS_SIZE, BASIS_SIZE));

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

		compute_quartics();
		fill_M_N();
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