#include <nlohmann/json.hpp>
#include <mrock/info.h>
#include <mrock/utility/OutputConvenience.hpp>
#include <mrock/utility/info_to_json.hpp>
#include <mrock/utility/FunctionTime.hpp>

#include <mrock/utility/Selfconsistency/BroydenSolver.hpp>

#include "DWave/T_C.hpp"
#include "../build_header/info.h"

using namespace DWave;
const std::string BASE_FOLDER = "../../data/dwave/";

int main(int argc, char** argv) {
    if (argc < 2) {
		std::cerr << "Invalid number of arguments: Use <path_to_executable> <configfile>" << std::endl;
		return -1;
	}
    mrock::utility::InputFileReader input(argv[1]);

	if (true) {
		Model model(input);
		auto solver = mrock::utility::Selfconsistency::make_broyden<l_float>(&model, &(model.Delta), 300);

		constexpr size_t BROYDEN_ITER = 700;
		constexpr double BROYDEN_EPS = 1e-8;
		solver.compute(false, BROYDEN_ITER, BROYDEN_EPS);

		std::cout << "########################################\n"
				  << "# Delta_max  = " << model.delta_max() << "\n"
				  << "# Delta_true = " << model.delta_true() << "\n"
				  << "########################################" << std::endl;

		const std::string output_folder = BASE_FOLDER + input.getString("output_folder") + "/" + model.to_folder();
		std::filesystem::create_directories(output_folder);
		std::cout << "Saving data to " << output_folder << std::endl;

		nlohmann::json comments = {
			{ "time", 				   	mrock::utility::time_stamp() },
			{ "V", 						model.dwave_coupling_in},
			{ "g", 					   	model.phonon_coupling_in },
		    { "E_F", 				   	model.fermi_energy },
			{ "omega_D", 			   	model.omega_debye },
		    { "N",              	   	model.N },
		    { "filling_at_zero_temp",	model.filling_at_zero_temp },
		};

		nlohmann::json jDelta = {
			{ "Delta",				model.Delta.as_vector(model.N * model.N) },
			{ "Delta_max",			model.delta_max() },
			{ "Delta_true",			model.delta_true() },
			{ "chemical_potential", model.chemical_potential }
		};

		jDelta.merge_patch(comments);
		mrock::utility::saveString(jDelta.dump(4), output_folder + "single_gap.json.gz");

		return 0;
	}

	// Compute T_C
	std::cout << "Starting T_C computation..." << std::endl;
	T_C tc(input);
	mrock::utility::member_function_time_ms(tc, &T_C::compute);

	nlohmann::json comments = {
		{ "time", 				   	mrock::utility::time_stamp() },
		{ "V", 						tc.model.dwave_coupling_in},
		{ "g", 					   	tc.model.phonon_coupling_in },
	    { "E_F", 				   	tc.model.fermi_energy },
		{ "omega_D", 			   	tc.model.omega_debye },
	    { "N",              	   	tc.model.N },
	    { "filling_at_zero_temp",	tc.model.filling_at_zero_temp },
	};
	
	nlohmann::json info_json = mrock::utility::generate_json<DWave::info>("dwave_");
	info_json.update(mrock::utility::generate_json<mrock::info>("mrock_"));

	const std::string output_folder = BASE_FOLDER + input.getString("output_folder") + "/" + tc.model.to_folder();
	std::filesystem::create_directories(output_folder);
	std::cout << "Saving data to " << output_folder << std::endl;
	mrock::utility::saveString(info_json.dump(4), output_folder + "metadata.json.gz");

	nlohmann::json jT_C = {
		{ "temperatures",			tc.temperatures },
		{ "delta_maxs", 			tc.delta_maxs },
		{ "delta_trues", 			tc.delta_trues },
		{ "chemical_potentials",	tc.chemical_potentials }
	};
	// All gaps are seperated because the data is rarely needed and takes a lot of time loading
	nlohmann::json jAllGaps = { { "finite_gaps",		tc.finite_gaps } };

	jT_C.merge_patch(comments);
	mrock::utility::saveString(jT_C.dump(4), output_folder + "T_C.json.gz");

	jAllGaps.merge_patch(comments);
	mrock::utility::saveString(jAllGaps.dump(4), output_folder + "all_gaps.json.gz");
	std::cout << "T_C computation finished." << std::endl;

    return 0;
}
