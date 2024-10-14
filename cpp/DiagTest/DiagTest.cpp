//#define EIGEN_USE_LAPACKE

#define LAPACK_COMPLEX_CUSTOM
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

#include <iostream>
#include <Eigen/Dense>
#include <omp.h>
#include <chrono>

using namespace Eigen;
using namespace std;

int main() {
#ifdef EIGEN_USE_LAPACKE
	cout << "Using LAPACK!" << endl;
#else
	cout << "Not using LAPACK!" << endl;
#endif
	const int N = 10000;
	srand(0); // Setting a seed
	MatrixXd B = MatrixXd::Random(N, N);
	MatrixXd A = 10 * B * B.transpose() / N;
	std::cout << "Matrix created!   " << A(0, 0) << std::endl;
	auto start = std::chrono::high_resolution_clock::now();
	SelfAdjointEigenSolver<MatrixXd> solver(A);
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	cout << "\nElapsed time for diagonalization: " << duration.count() << " milliseconds." << endl;

	cout << "Reconstruction error = "
		<< (A - solver.eigenvectors() * solver.eigenvalues().asDiagonal() * solver.eigenvectors().adjoint()).norm()
		<< endl;

	return 0;
}