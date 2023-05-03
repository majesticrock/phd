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

constexpr int siteIndex(int kx, int ky){
    return 2 * SYSTEM_SIZE * ky + kx;
}

constexpr long occupation(long identifier, int kx, int ky, bool spin){
    const int site_number = siteIndex(kx, ky);
    return (identifier >> (2 * site_number + ((spin) ? 1 : 0))) & 1;
}

constexpr long particleCreation(long identifier, int kx, int ky, bool spin){
    const int site_number = siteIndex(kx, ky);
    const short X = occupation(identifier, kx, ky, spin);
    if(X == 0){
        // Empty site
        identifier ^= 1UL << (2 * site_number + ((spin) ? 1 : 0));
        return identifier;
    }
    return -1;
}

constexpr long particleAnnihilation(long identifier, int kx, int ky, bool spin){
    const int site_number = siteIndex(kx, ky);
    const short X = occupation(identifier, kx, ky, spin);
    if(X == 1){
        // Filled site
        identifier ^= 1UL << (2 * site_number + ((spin) ? 1 : 0));
        return identifier;
    }
    return -1;
}

double unperturbed_energy(double kx, double ky){
    return -2*(cos(kx) + cos(ky));
}

void printSystem(long identifier){
    if(identifier < 0){
        std::cerr << "Invalid indentifier:   " << identifier << std::endl;
        return;
    }
    for (int i = -SYSTEM_SIZE; i < SYSTEM_SIZE; i++)
    {
        for (int j = -SYSTEM_SIZE; j < SYSTEM_SIZE; j++)
        {
            if(occupation(identifier, i + SYSTEM_SIZE, j + SYSTEM_SIZE, true)){
                std::cout << "U";
            } else {
                std::cout << "0";
            }
            if(occupation(identifier, i + SYSTEM_SIZE, j + SYSTEM_SIZE, false)){
                std::cout << "D";
            } else {
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
        for (int i = 0; i < 2*SYSTEM_SIZE; i++)
        {
            for (int j = 0; j < 2*SYSTEM_SIZE; j++)
            {
                long buffer = particleCreation(iden, i, j, true);
                if(buffer > -1){
                    creators[2 * siteIndex(i, j) + 1].coeffRef(buffer, iden) = 1;
                }
                buffer = particleCreation(iden, i, j, false);
                if(buffer > -1){
                    creators[2 * siteIndex(i, j)].coeffRef(buffer, iden) = 1;
                }
            }
        }
    }

    std::vector<SMatrix> annihilators;
    annihilators.reserve(OP_NUM);
    for(const auto& c : creators){
        annihilators.push_back(c.transpose());
    }

    const double U = -2.0 / (4 * SYSTEM_SIZE);
    const double V = 0.1 / (16 * 4 * SYSTEM_SIZE);
    const double mu = 0.5 * (U + V);

    SMatrix H(DIMENSION, DIMENSION);
    for (int i = 0; i < 2 * SYSTEM_SIZE; i++)
    {
        for (int j = 0; j < 2 * SYSTEM_SIZE; j++)
        {
            H += (unperturbed_energy((M_PI * i) / (2 * SYSTEM_SIZE) - M_PI, (M_PI * j) / (2 * SYSTEM_SIZE) - M_PI) - mu )
                * (creators[2 * siteIndex(i, j) + 1] * annihilators[2 * siteIndex(i, j) + 1]
                    + creators[2 * siteIndex(i, j)] * annihilators[2 * siteIndex(i, j)] );
        }
    }

    for (int i = 0; i < 2 * SYSTEM_SIZE; i++)
    {
        for (int j = 0; j < 2 * SYSTEM_SIZE; j++)
        {
            H += U * ( creators[2 * siteIndex(i, j) + 1] * creators[2 * siteIndex(i, j)] 
                * annihilators[2 * siteIndex(i, j)] * annihilators[2 * siteIndex(i, j) + 1]);
        }
    }
    
    for (int i = 0; i < 2 * SYSTEM_SIZE; i++)
    {
        for (int j = 0; j < 2 * SYSTEM_SIZE; j++)
        {
            for (int spin_l = 0; spin_l < 2; spin_l++)
            {
                for (int spin_r = 0; spin_r < 2; spin_r++)
                {
                    int i_r = (i + 1) % (2 * SYSTEM_SIZE);
                    int j_r = j;
                    H += V * (cos((M_PI * i) / (2 * SYSTEM_SIZE) - M_PI) + cos((M_PI * j) / (2 * SYSTEM_SIZE) - M_PI))
                        * ( creators[2 * siteIndex(i, j) + spin_l] * creators[2 * siteIndex(i_r, j_r) + spin_r] 
                        * annihilators[2 * siteIndex(i_r, j_r) + spin_r] * annihilators[2 * siteIndex(i, j) + spin_l]);
                    
                    if(i == 0){
                        i = 1;
                    } else {
                        i_r = (i - 1) % (2 * SYSTEM_SIZE);
                    }
                    j_r = j;
                    H += V * (cos((M_PI * i) / (2 * SYSTEM_SIZE) - M_PI) + cos((M_PI * j) / (2 * SYSTEM_SIZE) - M_PI))
                        * ( creators[2 * siteIndex(i, j) + spin_l] * creators[2 * siteIndex(i_r, j_r) + spin_r] 
                        * annihilators[2 * siteIndex(i_r, j_r) + spin_r] * annihilators[2 * siteIndex(i, j) + spin_l]);
                    
                    i_r = i;
                    j_r = (j + 1) % (2 * SYSTEM_SIZE);
                    H += V * (cos((M_PI * i) / (2 * SYSTEM_SIZE) - M_PI) + cos((M_PI * j) / (2 * SYSTEM_SIZE) - M_PI))
                        * ( creators[2 * siteIndex(i, j) + spin_l] * creators[2 * siteIndex(i_r, j_r) + spin_r] 
                        * annihilators[2 * siteIndex(i_r, j_r) + spin_r] * annihilators[2 * siteIndex(i, j) + spin_l]);
                    
                    i_r = i;
                    if(j == 0){
                        j = 1;
                    } else {
                        j_r = (j - 1) % (2 * SYSTEM_SIZE);
                    }
                    H += V * (cos((M_PI * i) / (2 * SYSTEM_SIZE) - M_PI) + cos((M_PI * j) / (2 * SYSTEM_SIZE) - M_PI))
                        * ( creators[2 * siteIndex(i, j) + spin_l] * creators[2 * siteIndex(i_r, j_r) + spin_r] 
                        * annihilators[2 * siteIndex(i_r, j_r) + spin_r] * annihilators[2 * siteIndex(i, j) + spin_l]);
                }
            }
        }
    }

    Eigen::SelfAdjointEigenSolver<SMatrix> solver(H);
    std::cout << solver.eigenvectors().col(0).dot(creators[2 * siteIndex(0, 0) + 1] * annihilators[2 * siteIndex(0, 0) + 1] * solver.eigenvectors().col(0)) << std::endl;
    return 0;
}

