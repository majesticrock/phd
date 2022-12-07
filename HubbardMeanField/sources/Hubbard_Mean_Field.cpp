#define _USE_MATH_DEFINES

#include <omp.h>
#include <mpi.h>

#include <string>
#include <iostream>
#include <filesystem>

#include "Utility/InputFileReader.hpp"
#include "Utility/OutputWriter.hpp"
#include "Hubbard/BasicHubbardModel.hpp"
#include "Hubbard/HubbardCDW.hpp"
#include "Hubbard/UsingBroyden.hpp"
#include "Hubbard/Constants.hpp"

int Hubbard::Constants::K_DISCRETIZATION = 100;

std::ostream& operator<<(std::ostream& os, Hubbard::Model::ModelParameters& mp) {
	os << mp.temperature << "\t" << mp.U << "\t" << mp.V;
	return os;
}

int main(int argc, char** argv)
{
	if (argc < 2) {
		std::cerr << "Invalid number of arguments: Use mpirun -n <threads> <path_to_executable> <configfile>" << std::endl;
		return -1;
	}
	// First call MPI_Init
	MPI_Init(&argc, &argv);

	// Get my rank and the number of ranks
	int rank, numberOfRanks;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfRanks);

	Utility::InputFileReader input(argv[1]);
	Hubbard::Constants::K_DISCRETIZATION = input.getInt("k_discretization");

//#define _DO_TEST
#ifdef _DO_TEST
	Hubbard::Model::ModelParameters mP(0, -1, 0.5, 0, 0, "", "");
	Hubbard::HubbardCDW model(mP);
	model.compute(true).print();
	std::cout << std::endl;
	Hubbard::UsingBroyden model2(mP);
	model2.compute(true).print();

	std::vector<std::vector<double>> energies;
	model2.getEnergies(energies, 1);
	Utility::saveData(energies, "../data/energies.txt");
	return MPI_Finalize();
#endif // _DO_TEST

	// Setup the parameters T, U, V
	std::vector<double> model_params = input.getDoubleList("model_parameters");

	// Setup the number of steps
	int GLOBAL_IT_STEPS = input.getInt("global_iterator_steps");
	int FIRST_IT_STEPS = GLOBAL_IT_STEPS / numberOfRanks;
	int SECOND_IT_STEPS = input.getInt("second_iterator_steps");

	double GLOBAL_IT_LIMS[2] = { 0, input.getDouble("global_iterator_upper_limit") };
	std::vector<std::string> option_list = { "T", "U", "V" };

	double SECOND_IT_MIN = 0, SECOND_IT_MAX = input.getDouble("second_iterator_upper_limit");
	double FIRST_IT_RANGE = 0;
	double FIRST_IT_MIN = 0, FIRST_IT_MAX = 0;
	for (int i = 0; i < option_list.size(); i++)
	{
		if (input.getString("global_iterator_type") == option_list[i]) {
			GLOBAL_IT_LIMS[0] = model_params[i];
			FIRST_IT_RANGE = (GLOBAL_IT_LIMS[1] - GLOBAL_IT_LIMS[0]) / numberOfRanks;
			FIRST_IT_MIN = GLOBAL_IT_LIMS[0] + rank * FIRST_IT_RANGE;
			FIRST_IT_MAX = FIRST_IT_MIN + FIRST_IT_RANGE;
			model_params[i] = FIRST_IT_MIN;
		}
		if (input.getString("second_iterator_type") == option_list[i]) {
			SECOND_IT_MIN = model_params[i];
		}
	}

	Hubbard::Model::ModelParameters modelParameters(model_params[0], model_params[1], model_params[2],
		(FIRST_IT_MAX - FIRST_IT_MIN) / FIRST_IT_STEPS, (SECOND_IT_MAX - SECOND_IT_MIN) / SECOND_IT_STEPS,
		input.getString("global_iterator_type"), input.getString("second_iterator_type"));

	typedef std::vector<double> data_vector;
	data_vector data_cdw(FIRST_IT_STEPS * SECOND_IT_STEPS);
	data_vector  data_sc(FIRST_IT_STEPS * SECOND_IT_STEPS);
	data_vector data_eta(FIRST_IT_STEPS * SECOND_IT_STEPS);

	for (int T = 0; T < FIRST_IT_STEPS; T++)
	{
#pragma omp parallel for num_threads(4)
		for (int U = 0; U < SECOND_IT_STEPS; U++)
		{
			Hubbard::Model::data_set ret;
			if (input.getBool("use_broyden")) {
				Hubbard::UsingBroyden model(modelParameters);
				ret = model.compute();
			}
			else {
				Hubbard::HubbardCDW model(modelParameters);
				ret = model.compute();
			}

			data_cdw[(T * SECOND_IT_STEPS) + U] = ret.delta_cdw;
			data_sc[(T * SECOND_IT_STEPS) + U] = ret.delta_sc;
			data_eta[(T * SECOND_IT_STEPS) + U] = ret.delta_eta;
			modelParameters.incrementSecondIterator();
		}

		modelParameters.printGlobal();
		std::cout << " done!" << std::endl;
		modelParameters.incrementGlobalIterator();
	}

	std::vector<double> recieve_cdw, recieve_sc, recieve_eta;
	if (rank == 0) {
		recieve_cdw.resize(GLOBAL_IT_STEPS * SECOND_IT_STEPS);
		recieve_sc.resize(GLOBAL_IT_STEPS * SECOND_IT_STEPS);
		recieve_eta.resize(GLOBAL_IT_STEPS * SECOND_IT_STEPS);
	}
	MPI_Gather(data_cdw.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, recieve_cdw.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(data_sc.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, recieve_sc.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(data_eta.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, recieve_eta.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		std::vector<std::string> comments;
		comments.push_back(input.getString("second_iterator_type") + "_MIN=" + std::to_string(SECOND_IT_MIN)
			+ "   " + input.getString("second_iterator_type") + "_MAX=" + std::to_string(SECOND_IT_MAX));

		comments.push_back(input.getString("global_iterator_type") + "_MIN=" + std::to_string(GLOBAL_IT_LIMS[0])
			+ "   " + input.getString("global_iterator_type") + "_MAX=" + std::to_string(GLOBAL_IT_LIMS[1]));

		std::string output_folder = input.getString("output_folder");
		std::filesystem::create_directories("../data/" + output_folder);

		Utility::saveData(recieve_cdw, SECOND_IT_STEPS, "../data/" + output_folder + "cdw.txt", comments);
		Utility::saveData(recieve_sc, SECOND_IT_STEPS, "../data/" + output_folder + "sc.txt", comments);
		Utility::saveData(recieve_eta, SECOND_IT_STEPS, "../data/" + output_folder + "eta.txt", comments);
	}

	return MPI_Finalize();
}