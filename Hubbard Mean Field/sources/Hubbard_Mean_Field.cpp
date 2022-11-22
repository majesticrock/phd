#define _USE_MATH_DEFINES

#include <iostream>
#include <omp.h>
#include <string>
#include "Utility/InputFileReader.hpp"
#include "Utility/OutputWriter.hpp"
#include "Hubbard/BasicHubbardModel.hpp"
#include "Hubbard/HubbardCDW.hpp"
#include "Hubbard/UsingBroyden.hpp"
#include "Hubbard/Constants.hpp"

int Hubbard::Constants::K_DISCRETIZATION = 100;

int main()
{
	omp_set_num_threads(16);

	Utility::InputFileReader input("params.config");
	Hubbard::Constants::K_DISCRETIZATION = input.getInt("k_discretization");

	const int T_STEPS = 64;
	const int U_STEPS = 128;

	typedef std::vector<double> data_vector;
	std::vector<data_vector> data_cdw(T_STEPS, data_vector(U_STEPS + 1));
	std::vector<data_vector> data_sc(T_STEPS, data_vector(U_STEPS + 1));
	std::vector<data_vector> data_eta(T_STEPS, data_vector(U_STEPS + 1));
	data_vector T_vals(T_STEPS);
	data_vector U_vals(U_STEPS + 1);

	const double T_MIN = 0, T_MAX = 1.2;
	const double U_MIN = -5, U_MAX = 0;

#define DO_TEST1
#ifdef DO_TEST
	double test_vals[] = { 1, -4.5 };
	Hubbard::HubbardCDW model(test_vals[0], test_vals[1], -0.5);
	model.compute(true).print();

	Hubbard::UsingBroyden model2(test_vals[0], test_vals[1], -0.5);
	model2.compute(true).print();
	return 0;
#endif // DO_TEST

#pragma omp parallel for schedule(dynamic)
	for (int T = 0; T < T_STEPS; T++)
	{
		double T_val = T_MIN + ((T_MAX - T_MIN) * T) / T_STEPS;
		T_vals[T] = T_val;
		for (int U = 0; U <= U_STEPS; U++)
		{
			double U_val = U_MIN + ((U_MAX - U_MIN) * U) / U_STEPS;
			U_vals[U] = U_val;
			//Hubbard::HubbardCDW model(T_val, U_val, -0.5);
			//Hubbard::HubbardCDW::data_set ret = model.compute();
			Hubbard::UsingBroyden model(T_val, U_val, -0.5);
			Hubbard::UsingBroyden::data_set ret = model.compute();

			data_cdw[T][U] = ret.delta_cdw;
			data_sc[T][U] = ret.delta_sc;
			data_eta[T][U] = ret.delta_eta;
		}

		std::cout << "T=" << T_val << " done!" << std::endl;
	}
	std::vector<std::string> comments;
	comments.push_back("U_min=" + std::to_string(U_MIN) + "   U_max=" + std::to_string(U_MAX));
	comments.push_back("T_min=" + std::to_string(T_MIN) + "   T_max=" + std::to_string(T_MAX));

	Utility::saveData(data_cdw, "../data/basic_hubbard_cdw.txt", comments);
	Utility::saveData(data_sc, "../data/basic_hubbard_sc.txt", comments);
	Utility::saveData(data_eta, "../data/basic_hubbard_eta.txt", comments);

	return 0;
}