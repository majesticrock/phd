/* This file serves as a testing environment for the various
* algorithms provided in this library
*/
#define _USE_MATH_DEFINES

#include "Numerics/MidpointRule.hpp"
#include "Numerics/TrapezoidalRule.hpp"

#include <cmath>
#include <iostream>

using namespace Utility;
using namespace Utility::Numerics;

int main(int argc, char** argv) {
	auto test_func = [](double t, double f) {
		return cos(M_PI * f * t);
		};
	auto test_func_2 = [](double t) {
		return cos(M_PI * t);
		};

	constexpr int step_num = 1000;

	std::cout << "Test 1: f=1:\n";
	std::cout << Integration::midpoint_rule(test_func, 0., 5., step_num, 1) << "\t"
		<< Integration::trapezoidal_rule(test_func, 0., 5., step_num, 1) << std::endl;

	std::cout << "Test 1: f=0.1:\n";
	std::cout << Integration::midpoint_rule(test_func, 0., 2., step_num, 0.1) << "\t"
		<< Integration::trapezoidal_rule(test_func, 0., 2., step_num, 0.1) << std::endl;

	std::cout << "Test 2: f=0.1:\n";
	std::cout << Integration::midpoint_rule(test_func_2, 0., 0.3, step_num) << "\t"
		<< Integration::trapezoidal_rule(test_func_2, 0., 0.3, step_num) << std::endl;

	return 0;
}