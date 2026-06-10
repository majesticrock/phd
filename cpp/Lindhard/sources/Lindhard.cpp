#include <mrock/utility/OutputConvenience.hpp>
#include <chrono>
#include <iostream>
#include <omp.h>
#include <nlohmann/json.hpp>
#include "helper_functions.hpp"

std::vector<cdouble> lindhard(const kvec_t& q, 
                const std::vector<cdouble>& omegas,
                const triple_array<kvec_t>& cosines,
                const triple_array<kvec_t>& sines,
                const triple_array<double>& xis,
                const triple_array<double>& fermis)
{
    const double qnorm_sqr = q[0]*q[0] + q[1]*q[1] + q[2]*q[2];
    std::vector<cdouble> chis(omegas.size());

    for (size_t o=0U; o<omegas.size(); ++o) {
        if (omegas[o].real() < omegas[o].imag() && qnorm_sqr < omegas[o].imag()*omegas[o].imag()) {
            // help
            for (size_t x=0U; x<N; ++x) {
                for (size_t y=0U; y<N; ++y) {
                    for (size_t z=0U; z<N; ++z) {
                        chis[o] += derivative_fermi(xis[x][y][z]);
                    }
                }
            }
        }
        else {
            if (qnorm_sqr < omegas[o].imag()*omegas[o].imag())
            {
                for (size_t x=0U; x<N; ++x) {
                    for (size_t y=0U; y<N; ++y) {
                        for (size_t z=0U; z<N; ++z) {
                            chis[o] += fermis[x][y][z] * (shifted_xi(q, cosines[x][y][z], sines[x][y][z]) - xis[x][y][z]);
                        }
                    }
                }
                chis[o] /= omegas[o]*omegas[o];
            }
            else {
                double xi_q;
                for (size_t x=0U; x<N; ++x) {
                    for (size_t y=0U; y<N; ++y) {
                        for (size_t z=0U; z<N; ++z) {
                            xi_q = shifted_xi(q, cosines[x][y][z], sines[x][y][z]); 
                            chis[o] += (fermis[x][y][z] - fermi(xi_q)) / (xis[x][y][z] - xi_q + omegas[o]);
                        }
                    }
                }
            }
        }
        chis[o] *= -2. / (N*N*N);
    }
    return chis;
}

int main(int argc, char** argv) {
    using clock = std::chrono::high_resolution_clock;
    clock::time_point begin = clock::now();

    std::vector<cdouble> omegas = {
         {0.0, 1e-5},
         {0.1, 5e-3},
         {0.2, 5e-3},
         {0.5, 5e-3},
    };

    const triple_array<kvec_t> kvecs = generate_kvecs();
    const triple_array<kvec_t> cosines = compute_cosines(kvecs);;
    const triple_array<kvec_t> sines = compute_sines(kvecs);;
    const triple_array<double> xis = compute_xis(cosines);
    const triple_array<double> fermis = compute_fermis(xis);

    std::array<kvec_t, 4> high_symmetry_points{}; // Gamma, X, M, R
    // Gamma
    high_symmetry_points[0] = {0., 0., 0.};
    // X
    high_symmetry_points[1] = multiply(0.5, basis_vectors[0]);
    // M
    high_symmetry_points[2] = multiply(0.5, add(basis_vectors[0], basis_vectors[1]));
    // R
    high_symmetry_points[3] = multiply(0.5, add(basis_vectors[0], add(basis_vectors[1], basis_vectors[2])));

    std::array<std::vector<std::array<double, n_segment>>, 4> chi_paths;
    chi_paths.fill(std::vector<std::array<double, n_segment>>(omegas.size()));

    for (size_t i=0U; i<4U; ++i) {
        const auto& qvecs = connector_line(high_symmetry_points[i], high_symmetry_points[(i+1)%4]);
#pragma omp parallel for
        for (int q=0; q<n_segment; ++q) {
            std::vector<cdouble> result = lindhard(qvecs[q], omegas, cosines, sines, xis, fermis);
            for (size_t r=0U; r<result.size(); ++r) {
                chi_paths[i][r][q] = result[r].real();
            }
        }
    }

    /*
    * output  data
    */
    std::vector<double> real_omegas(omegas.size());
    for(size_t o=0U;o<omegas.size();++o){
        real_omegas[o] = omegas[o].real();
    }
	nlohmann::json jChi = {
		{ "chi_Gamma_X",	chi_paths[0] },
	    { "chi_X_M",    	chi_paths[1] },
        { "chi_M_R",    	chi_paths[2] },
        { "chi_R_Gamma",	chi_paths[3] },
        { "omegas",         real_omegas  }
	};
	mrock::utility::saveString(jChi.dump(4), "chi.json.gz");
	std::cout << "Data saved!" << std::endl;

    clock::time_point end = clock::now();
	std::cout << "Runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    return 0;
}