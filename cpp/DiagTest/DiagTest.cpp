#include <iostream>
#include <armadillo>
#include <Eigen/Dense>
#include <chrono>

int main() {
    const int n = 4000; // Matrix size

    // Generate a symmetric matrix using Eigen
    Eigen::MatrixXd x = Eigen::MatrixXd::Random(n, n);
    Eigen::MatrixXd matrix = (x + x.transpose()) / 2; // Make sure it's symmetric

    // ----------------- Armadillo -----------------
    // Start measuring time for Armadillo
    auto start_arma = std::chrono::high_resolution_clock::now();

    // Convert Eigen matrix to Armadillo matrix (sharing memory)
    arma::mat A = arma::mat(matrix.data(), matrix.rows(), matrix.cols(), false, true);

    // Perform eigenvalue decomposition with Armadillo
    arma::vec eigval_arma;
    arma::mat eigvec_arma;
    arma::eig_sym(eigval_arma, eigvec_arma, A, "std");  // eig_sym computes for symmetric matrices

    // Reconstruct matrix from Armadillo eigenvalues and eigenvectors
    arma::mat reconstructed_arma = eigvec_arma * arma::diagmat(eigval_arma) * eigvec_arma.t();

    // End measuring time for Armadillo
    auto end_arma = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_arma = end_arma - start_arma;
    std::cout << "Armadillo decomposition took " << elapsed_arma.count() << " seconds." << std::endl;

    // Compute norm of the difference (original - reconstructed)
    Eigen::MatrixXd reconstructed_arma_eigen = Eigen::Map<Eigen::MatrixXd>(reconstructed_arma.memptr(), n, n);
    double arma_diff_norm = (matrix - reconstructed_arma_eigen).norm();
    std::cout << "Norm of difference (original - reconstructed) using Armadillo: " << arma_diff_norm << std::endl;

    // ----------------- Eigen -----------------
    // Start measuring time for Eigen
    auto start_eigen = std::chrono::high_resolution_clock::now();

    // Perform eigenvalue decomposition using Eigen
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);
    Eigen::VectorXd eigenvalues_eigen = solver.eigenvalues();
    Eigen::MatrixXd eigenvectors_eigen = solver.eigenvectors();

    // Reconstruct matrix from Eigen eigenvalues and eigenvectors
    Eigen::MatrixXd reconstructed_eigen = eigenvectors_eigen * eigenvalues_eigen.asDiagonal() * eigenvectors_eigen.transpose();

    // End measuring time for Eigen
    auto end_eigen = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_eigen = end_eigen - start_eigen;
    std::cout << "Eigen decomposition took " << elapsed_eigen.count() << " seconds." << std::endl;

    // Compute norm of the difference (original - reconstructed)
    double eigen_diff_norm = (matrix - reconstructed_eigen).norm();
    std::cout << "Norm of difference (original - reconstructed) using Eigen: " << eigen_diff_norm << std::endl;

    return 0;
}
