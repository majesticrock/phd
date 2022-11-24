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

int main(int argc, char** argv)
{
	// First call MPI_Init
	MPI_Init(&argc, &argv);

	// Get my rank and the number of ranks
	int rank, numberOfRanks;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfRanks);

	Utility::InputFileReader input("params.config");
	Hubbard::Constants::K_DISCRETIZATION = input.getInt("k_discretization");

	const int GLOBAL_T_STEPS = 160;
	int T_STEPS = GLOBAL_T_STEPS / numberOfRanks;
	int U_STEPS = 160;
	// T_MIN, T_MAX for all ranks
	const double GLOBAL_T_LIMS[4] = {-2, 2};
	// Limits for the current rank
	double T_RANK_RANGE = (GLOBAL_T_LIMS[1] - GLOBAL_T_LIMS[0]) / numberOfRanks;
	double T_MIN = GLOBAL_T_LIMS[0] + rank * T_RANK_RANGE, T_MAX = T_MIN + T_RANK_RANGE;
	double U_MIN = -5, U_MAX = 0;

	typedef std::vector<double> data_vector;
	data_vector data_cdw(T_STEPS * (U_STEPS + 1));
	data_vector  data_sc(T_STEPS * (U_STEPS + 1));
	data_vector data_eta(T_STEPS * (U_STEPS + 1));

//#define _DO_TEST
#ifdef _DO_TEST
	double test_vals[] = { 0, 0, -0.025 };
	Hubbard::HubbardCDW model(test_vals[0], test_vals[1], test_vals[2]);
	model.compute(true).print();
	std::cout << std::endl;
	Hubbard::UsingBroyden model2(test_vals[0], test_vals[1], test_vals[2]);
	model2.compute(true).print();
	return MPI_Finalize();
#endif // _DO_TEST

	for (int T = 0; T < T_STEPS; T++)
	{
		double T_val = T_MIN + ((T_MAX - T_MIN) * T) / T_STEPS;
		//#pragma omp parallel for schedule(dynamic)
		for (int U = 0; U <= U_STEPS; U++)
		{
			double U_val = U_MIN + ((U_MAX - U_MIN) * U) / U_STEPS;
			//Hubbard::BasicHubbardModel model(T_val, U_val);
			//Hubbard::BasicHubbardModel::data_set ret = model.compute();
			Hubbard::UsingBroyden model(0, U_val, T_val);
			Hubbard::UsingBroyden::data_set ret = model.compute();

			data_cdw[(T * (U_STEPS + 1)) + U] = ret.delta_cdw;
			data_sc[(T * (U_STEPS + 1)) + U] = ret.delta_sc;
			data_eta[(T * (U_STEPS + 1)) + U] = ret.delta_eta;
		}

		std::cout << "T=" << T_val << " done!" << std::endl;
	}

	std::vector<double> recieve_cdw, recieve_sc, recieve_eta;
	if(rank == 0){
		recieve_cdw.resize(GLOBAL_T_STEPS * (U_STEPS + 1));
		recieve_sc.resize(GLOBAL_T_STEPS * (U_STEPS + 1));
		recieve_eta.resize(GLOBAL_T_STEPS * (U_STEPS + 1));
	}
	MPI_Gather(data_cdw.data(), T_STEPS * (U_STEPS + 1), MPI_DOUBLE, recieve_cdw.data(), T_STEPS * (U_STEPS + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(data_sc.data(), T_STEPS * (U_STEPS + 1), MPI_DOUBLE, recieve_sc.data(), T_STEPS * (U_STEPS + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(data_eta.data(), T_STEPS * (U_STEPS + 1), MPI_DOUBLE, recieve_eta.data(), T_STEPS * (U_STEPS + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(rank == 0){
		std::vector<std::string> comments;
		comments.push_back("U_min=" + std::to_string(U_MIN) + "   U_max=" + std::to_string(U_MAX));
		comments.push_back("T_min=" + std::to_string(GLOBAL_T_LIMS[0]) + "   T_max=" + std::to_string(GLOBAL_T_LIMS[1]));

		std::string output_folder = input.getString("output_folder");
		std::filesystem::create_directories("../data/" + output_folder);
		
		Utility::saveData(recieve_cdw, U_STEPS + 1, "../data/" + output_folder + "cdw.txt", comments);
		Utility::saveData(recieve_sc, U_STEPS + 1 , "../data/" + output_folder + "sc.txt", comments);
		Utility::saveData(recieve_eta, U_STEPS + 1, "../data/" + output_folder + "eta.txt", comments);
	}

	return MPI_Finalize();
}