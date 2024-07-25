#include <omp.h>
#include <cmath>
#include <vector>
#include <chrono>
#include <iostream>
#include <mpi.h>

int main(int argc, char** argv){
	// First call MPI_Init
	MPI_Init(&argc, &argv);

	std::chrono::steady_clock::time_point test_b = std::chrono::steady_clock::now();
	std::chrono::steady_clock::time_point test_e;

	const std::vector<int> ranges = {
		100, 200, 500, 1000,
		2000, 5000, 10000,
		20000, 50000, 100000
	};

	for (auto range : ranges) {
		test_b = std::chrono::steady_clock::now();

		auto c = [](double k) {
			return cos(k);
		};

		std::vector<double> v(range);
#pragma omp parallel for
		for (size_t i = 0U; i < range; ++i)
		{
			for (size_t j = 0U; j < range; ++j)
			{
				v[i] += c(j);
			}
		}
		test_e = std::chrono::steady_clock::now();
		std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(test_e - test_b).count() << "[ms]" << std::endl;
	}

	return MPI_Finalize();
}