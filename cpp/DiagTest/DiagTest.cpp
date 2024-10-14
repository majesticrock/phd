#include <iostream>
#include <Eigen/Dense>
#include <omp.h>
#include <chrono>

using namespace Eigen;
using namespace std;

void householderTridiagonalization(const MatrixXd& input, MatrixXd& Q, VectorXd& d, VectorXd& e) {
    MatrixXd A = input;
    int n = A.rows();

    // Initialize Q to the identity matrix
    Q = MatrixXd::Identity(n, n);

    // Prepare diagonal and off-diagonal storage
    d.resize(n);
    e.resize(n - 1);

    for (int k = 0; k < n - 2; ++k) {
        // Extract the (k+1)-th column from the submatrix A(k+1:n, k)
        VectorXd x = A.block(k + 1, k, n - k - 1, 1);

        // Compute the Householder vector v
        VectorXd v = x;
        v(0) += (x(0) >= 0 ? 1.0 : -1.0) * x.norm();

        // Normalize the Householder vector
        v.normalize();

        // Apply the transformation A = (I - 2*v*v^T) * A * (I - 2*v*v^T)
        MatrixXd P = MatrixXd::Identity(n - k - 1, n - k - 1) - 2.0 * v * v.transpose();
            A.block(k + 1, k + 1, n - k - 1, n - k - 1) = P * A.block(k + 1, k + 1, n - k - 1, n - k - 1) * P;
            A.block(k + 1, k, n - k - 1, 1) = P * A.block(k + 1, k, n - k - 1, 1);
            A.block(k, k + 1, 1, n - k - 1) = A.block(k + 1, k, n - k - 1, 1).transpose();

        // Update the matrix Q
        MatrixXd Qk = MatrixXd::Identity(n, n);
        Qk.block(k + 1, k + 1, n - k - 1, n - k - 1) -= 2.0 * v * v.transpose();

            Q = Q * Qk;

        // Store the diagonal and off-diagonal elements
        d(k) = A(k, k);
        e(k) = A(k + 1, k);
    }

    // Store the last two diagonal elements
    d(n - 2) = A(n - 2, n - 2);
    d(n - 1) = A(n - 1, n - 1);
    e(n - 2) = A(n - 1, n - 2);
}

int main() {
    // Example: a random Hermitian matrix
    MatrixXd B = MatrixXd::Random(10000, 10000);
    MatrixXd A = B + B.transpose();  // Make it symmetric (Hermitian in the real case)

    auto start = std::chrono::high_resolution_clock::now();
    Tridiagonalization<MatrixXd> eig_t(A);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "\nElapsed time for tridiagonalization: " << duration.count() << " milliseconds." << endl;

    return 0;
}
