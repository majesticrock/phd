#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <cmath>
#include "../../HubbardMeanField/sources/Utility/OutputConvenience.hpp"

const double CUT_OFF = 1e-6;

void loadData(const std::string& filename, std::vector<double>& reciever) {
	std::string line;
	std::ifstream reader(filename);
	if (!reader.is_open()) throw std::runtime_error("Could not open file " + filename);
	while (std::getline(reader, line)) {
		if (line.empty()) continue;
		if (line.front() == '#') continue;

		std::istringstream iss(line);
		std::string data_point;
		while (std::getline(iss, data_point, ' ')) {
			reciever.emplace_back(std::stod(data_point));
		}
	}
	reader.close();
}

int main()
{
	const std::string filepath = "../../data/T0/U_modes/";
	const std::string filename = "-2.00";
	const std::string systems[] = { "10/" , "12/", "14/", "16/", "18/", "20/", "24/", "30/", "36/", "40/", "46/", "50/", "60/" };

	std::chrono::time_point begin = std::chrono::steady_clock::now();
	std::chrono::time_point end = std::chrono::steady_clock::now();
	for (const auto& system : systems) {
		begin = std::chrono::steady_clock::now();
		std::vector<double> one_particle_data;
		loadData(filepath + system + filename + "_one_particle.txt", one_particle_data);

		std::set<double> sum_of_two_particles;
		for (const double& eps_1 : one_particle_data) {
			for (const double& eps_2 : one_particle_data) {
				sum_of_two_particles.insert(eps_1 + eps_2);
			}
		}

		std::vector<double> two_particle_data;
		loadData(filepath + system + filename + ".txt", two_particle_data);

		std::vector<double> special_modes;
		for (const double& eps : two_particle_data) {
			if (std::abs(eps) < 1e-2) continue;
			auto it = sum_of_two_particles.lower_bound(1. / eps);
			double val_left = std::abs((*it) - 1. / eps);
			if (std::abs(val_left) < CUT_OFF) continue;

			if (++it == sum_of_two_particles.end()) {
				special_modes.push_back(val_left);
				continue;
			}
			double val_right = std::abs((*it) - 1. / eps);
			special_modes.push_back((std::abs(val_left) > std::abs(val_right)) ? val_right : val_left);
		}

		std::sort(special_modes.begin(), special_modes.end());

		end = std::chrono::steady_clock::now();
		std::cout << "Time: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

		Utility::saveData(special_modes, filepath + system + filename + "_special.txt");
		std::cout << "Done with system " << system << std::endl;
	}

	return 0;
}