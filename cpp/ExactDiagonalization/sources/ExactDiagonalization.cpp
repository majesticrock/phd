#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>

constexpr int SYSTEM_SIZE = 1;
constexpr int OP_NUM = 8 * SYSTEM_SIZE * SYSTEM_SIZE;
constexpr long DIMENSION = 1 << OP_NUM;

typedef Eigen::SparseMatrix<double> SMatrix;
typedef Eigen::MatrixXd DMatrix;

constexpr int siteIndex(int kx, int ky) {
	return 2 * SYSTEM_SIZE * ky + kx;
}

constexpr long occupation(long identifier, int kx, int ky, bool spin) {
	const int site_number = siteIndex(kx, ky);
	return (identifier >> (2 * site_number + ((spin) ? 1 : 0))) & 1;
}

constexpr long particleCreation(long identifier, int kx, int ky, bool spin) {
	const int site_number = siteIndex(kx, ky);
	const short X = occupation(identifier, kx, ky, spin);
	if (X == 0) {
		// Empty site
		identifier ^= 1UL << (2 * site_number + ((spin) ? 1 : 0));
		return identifier;
	}
	return -1;
}

constexpr long particleAnnihilation(long identifier, int kx, int ky, bool spin) {
	const int site_number = siteIndex(kx, ky);
	const short X = occupation(identifier, kx, ky, spin);
	if (X == 1) {
		// Filled site
		identifier ^= 1UL << (2 * site_number + ((spin) ? 1 : 0));
		return identifier;
	}
	return -1;
}

double cos_x_cos_y(int kx, int ky) {
	return cos((M_PI * kx) / SYSTEM_SIZE) + cos((M_PI * ky) / SYSTEM_SIZE);
}

double unperturbed_energy(int kx, int ky) {
	return -2 * cos_x_cos_y(kx, ky);
}

void printSystem(long identifier) {
	if (identifier < 0) {
		std::cerr << "Invalid indentifier:   " << identifier << std::endl;
		return;
	}
	for (int i = -SYSTEM_SIZE; i < SYSTEM_SIZE; i++)
	{
		for (int j = -SYSTEM_SIZE; j < SYSTEM_SIZE; j++)
		{
			if (occupation(identifier, i + SYSTEM_SIZE, j + SYSTEM_SIZE, true)) {
				std::cout << "U";
			}
			else {
				std::cout << "0";
			}
			if (occupation(identifier, i + SYSTEM_SIZE, j + SYSTEM_SIZE, false)) {
				std::cout << "D";
			}
			else {
				std::cout << "0";
			}
			std::cout << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

int main()
{
	std::vector<SMatrix> creators(OP_NUM, SMatrix(DIMENSION, DIMENSION));

	for (long iden = 0; iden < DIMENSION; iden++)
	{
		for (int i = 0; i < 2 * SYSTEM_SIZE; i++)
		{
			for (int j = 0; j < 2 * SYSTEM_SIZE; j++)
			{
				long buffer = particleCreation(iden, i, j, true);
				if (buffer > -1) {
					creators[2 * siteIndex(i, j) + 1].coeffRef(buffer, iden) = 1;
				}
				buffer = particleCreation(iden, i, j, false);
				if (buffer > -1) {
					creators[2 * siteIndex(i, j)].coeffRef(buffer, iden) = 1;
				}
			}
		}
	}

	std::vector<SMatrix> annihilators;
	annihilators.reserve(OP_NUM);
	for (const auto& c : creators) {
		annihilators.push_back(c.transpose());
	}

	constexpr double U = -2.0 / (4 * SYSTEM_SIZE * SYSTEM_SIZE);
	constexpr double V = 0 / (16 * 4 * SYSTEM_SIZE * SYSTEM_SIZE);
	constexpr double mu = 0.5 * (U + V) * 4 * SYSTEM_SIZE * SYSTEM_SIZE;

	SMatrix H(DIMENSION, DIMENSION);
	for (int i = 0; i < 2 * SYSTEM_SIZE; i++)
	{
		for (int j = 0; j < 2 * SYSTEM_SIZE; j++)
		{
			H += (unperturbed_energy(i - 1, j - 1) - mu)
				* (creators[2 * siteIndex(i, j) + 1] * annihilators[2 * siteIndex(i, j) + 1]
					+ creators[2 * siteIndex(i, j)] * annihilators[2 * siteIndex(i, j)]);
		}
	}

	auto fix_index = [](int index) -> int {
		if (index < 0) index += 2 * SYSTEM_SIZE;
		index %= SYSTEM_SIZE;
		return index;
	};

	for (int ki = 0; ki < 2 * SYSTEM_SIZE; ki++)
	{
		for (int kj = 0; kj < 2 * SYSTEM_SIZE; kj++)
		{
			for (int qi = 0; qi < 2 * SYSTEM_SIZE; qi++)
			{
				for (int qj = 0; qj < 2 * SYSTEM_SIZE; qj++)
				{
					for (int pi = 0; pi < 2 * SYSTEM_SIZE; pi++)
					{
						for (int pj = 0; pj < 2 * SYSTEM_SIZE; pj++)
						{
							H += U * (creators[2 * siteIndex(ki, kj) + 1] * creators[2 * siteIndex(pi, pj)]
								* annihilators[2 * siteIndex(fix_index(pi + qi), fix_index(pj + qi))]
								* annihilators[2 * siteIndex(fix_index(ki - qi), fix_index(kj - qi)) + 1]);
						}
					}
				}
			}
		}
	}

	Eigen::SelfAdjointEigenSolver<SMatrix> solver(H);
	int smallest = 0;
	for (size_t i = 1; i < solver.eigenvalues().size(); i++)
	{
		if (solver.eigenvalues()(i) < solver.eigenvalues()(smallest)) {
			smallest = i;
		}
	}
	SMatrix total_occupation(DIMENSION, DIMENSION);
	for (int i = 0; i < 2 * SYSTEM_SIZE; i++)
	{
		for (int j = 0; j < 2 * SYSTEM_SIZE; j++)
		{
			total_occupation += creators[2 * siteIndex(i, j) + 1] * annihilators[2 * siteIndex(i, j) + 1];
		}
	}

	std::cout << solver.eigenvectors().col(smallest).dot(
		total_occupation * solver.eigenvectors().col(smallest)) << std::endl;
	//std::cout << solver.eigenvectors().col(smallest) << std::endl;
	//std::cout << solver.eigenvalues() << std::endl;
	for (size_t i = 0; i < DIMENSION; i++)
	{
		if (std::abs(solver.eigenvectors()(i, smallest)) > 1e-8) {
			//std::cout << i << ":\n";
			//printSystem(i);
		}
	}

	return 0;
}