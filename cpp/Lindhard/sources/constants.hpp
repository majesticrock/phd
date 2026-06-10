#include <vector>
#include <array>
#include <complex>
#include <cmath>
#include <numbers>

typedef std::array<double,3> kvec_t;
typedef std::complex<double> cdouble;

// Constants
constexpr double pi = std::numbers::pi;
constexpr double a = 1.0;
constexpr double beta = 30.;
constexpr double W = 1.0;
constexpr double E_F = -0.5;
constexpr int N = 200;
constexpr int n_segment = 100;

template<class T>
using N_array = std::array<T, N>; 
template<class T>
using triple_array = N_array<N_array<N_array<T>>>;

constexpr std::array<kvec_t, 3> basis_vectors = {
    kvec_t{0.,       2*pi/a, 2*pi/a},
    kvec_t{2*pi/a, 0.,       2*pi/a},
    kvec_t{2*pi/a, 2*pi/a, 0.      }
};

constexpr kvec_t& addInPlace(kvec_t& k1, const kvec_t& k2) 
{
    for(size_t i=0U; i<k1.size(); ++i){
        k1[i] += k2[i];
    }
    return k1;
}
constexpr kvec_t& multiplyInPlace(const double scalar, kvec_t& kvec) 
{
    for (auto& k : kvec) {
        k *= scalar;
    }
    return kvec;
}

constexpr kvec_t multiply(const double scalar, kvec_t kvec) 
{
    return multiplyInPlace(scalar, kvec);
}
constexpr kvec_t add(kvec_t k1, const kvec_t& k2) 
{
    return addInPlace(k1, k2);
}

constexpr double coeff(int i) {
    return static_cast<double>(i - N/2) / static_cast<double>(N);
}
constexpr N_array<double> generate_coeffs() {
    N_array<double> coeffs;
    for (size_t i=0U; i<N; ++i) {
        coeffs[i] = coeff(i);
    }
    return coeffs;
}

constexpr triple_array<kvec_t> generate_kvecs () 
{
    constexpr N_array<double> coeffs = generate_coeffs();
    triple_array<kvec_t>kvecs{};

    for (size_t x=0U; x<N; ++x) {
        for (size_t y=0U; y<N; ++y) {
            for (size_t z=0U; z<N; ++z) {
                addInPlace(kvecs[x][y][z], multiply(coeffs[x], basis_vectors[0]));
                addInPlace(kvecs[x][y][z], multiply(coeffs[y], basis_vectors[1]));
                addInPlace(kvecs[x][y][z], multiply(coeffs[z], basis_vectors[2]));
            }
        }
    }
    return kvecs;
}
