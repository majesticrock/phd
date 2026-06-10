#include "constants.hpp"

inline triple_array<kvec_t> compute_cosines(const triple_array<kvec_t>& kvecs) 
{
    triple_array<kvec_t> cosines{};
    for (size_t x=0U; x<N; ++x) {
        for (size_t y=0U; y<N; ++y) {
            for (size_t z=0U; z<N; ++z) {
                for (size_t k=0; k<3; ++k) {
                    cosines[x][y][z][k] = std::cos(0.5 * a * kvecs[x][y][z][k]);
                }
            }
        }
    }
    return cosines;
}

inline triple_array<kvec_t> compute_sines(const triple_array<kvec_t>& kvecs) 
{
    triple_array<kvec_t> sines{};
    for (size_t x=0U; x<N; ++x) {
        for (size_t y=0U; y<N; ++y) {
            for (size_t z=0U; z<N; ++z) {
                for (size_t k=0; k<3; ++k) {
                    sines[x][y][z][k] = std::sin(0.5 * a * kvecs[x][y][z][k]);
                }
            }
        }
    }
    return sines;
}

inline triple_array<double> compute_xis(const triple_array<kvec_t>& cosines) 
{
    triple_array<double> xis{};
    for (size_t x=0U; x<N; ++x) {
        for (size_t y=0U; y<N; ++y) {
            for (size_t z=0U; z<N; ++z) {
                xis[x][y][z] = W;
                for (size_t k=0; k<3; ++k) {
                    xis[x][y][z] *= -cosines[x][y][z][k];
                }
                xis[x][y][z] -= E_F;
            }
        }
    }
    return xis;
}

inline double shifted_xi(const kvec_t& q, 
                         const kvec_t& cosines, 
                         const kvec_t& sines) 
{
    const kvec_t cosq = { std::cos(0.5*a*q[0]), std::cos(0.5*a*q[1]), std::cos(0.5*a*q[2]) };
    const kvec_t sinq = { std::sin(0.5*a*q[0]), std::sin(0.5*a*q[1]), std::sin(0.5*a*q[2]) };

    double xi_q{W};
    for (size_t k=0; k<3; ++k) {
        xi_q *= -(cosq[k] * cosines[k] - sinq[k] * sines[k]);
    }
    xi_q -= E_F;

    return xi_q;
}

inline triple_array<double> compute_fermis(const triple_array<double>& xis) 
{
    triple_array<double> fermis{};
    for (size_t x=0U; x<N; ++x) {
        for (size_t y=0U; y<N; ++y) {
            for (size_t z=0U; z<N; ++z) {
                fermis[x][y][z] = 1.0 / (std::exp(beta * xis[x][y][z]) + 1.0);
            }
        }
    }
    return fermis;
}

inline double fermi(const double x) {
    return 1.0 / (1.0 + std::exp(beta * x));
}
inline double derivative_fermi(const double x) {
    return -0.5 * beta / (1. + std::cosh(beta * x));
}

inline std::array<kvec_t, n_segment> connector_line(const kvec_t& start, const kvec_t& end) {
    const kvec_t connect = {end[0] - start[0], end[1] - start[1], end[2] - start[2]};
    std::array<kvec_t, n_segment> connector_line{};

    for(int i=0; i<n_segment; ++i) {
        connector_line[i] = add(start, multiply(static_cast<double>(i)/static_cast<double>(n_segment), connect));
    }
    return connector_line;
}