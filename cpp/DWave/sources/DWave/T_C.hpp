#pragma once

#include "GlobalDefinitions.hpp"
#include "Model.hpp"

#include <mrock/utility/InputFileReader.hpp>
#include <vector>

namespace DWave {
    struct T_C {
        Model model;
        std::vector<l_float> temperatures;
        std::vector<std::vector<l_float>> finite_gaps;
        std::vector<l_float> delta_maxs;
        std::vector<l_float> delta_trues;
        std::vector<l_float> chemical_potentials;

        T_C(mrock::utility::InputFileReader& input);

        void compute();
    };
}