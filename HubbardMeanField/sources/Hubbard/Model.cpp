#include "Model.hpp"
#include <iostream>

namespace Hubbard {
	void Model::data_set::print() const {
		std::cout << delta_cdw << "\t" << delta_sc << "\t" << delta_eta
			<< "\t" << sqrt(delta_cdw * delta_cdw + delta_sc * delta_sc + delta_eta * delta_eta) << std::endl;
	}

	Model::Model(double _temperature, double _U)
		: temperature(_temperature), U(_U)
	{
	}
	Model::Model(ModelParameters& _params)
		: temperature(_params.temperature), U(_params.U)
	{
	}
	void Model::ModelParameters::incrementer(std::string& s, const double step)
	{
		if (s == "T") {
			temperature += step;
		}
		else if (s == "U") {
			U += step;
		}
		else if (s == "V") {
			V += step;
		}
	}
	Model::ModelParameters::ModelParameters(double _temperature, double _U, double _V, double _global_step, double _second_step,
		std::string _global_iterator_type, std::string _second_iterator_type)
		: global_iterator_type(_global_iterator_type), second_iterator_type(_second_iterator_type),
		global_step(_global_step), second_step(_second_step), temperature(_temperature), U(_U), V(_V)
	{
		if (second_iterator_type == "T") {
			second_it_min = temperature;
		}
		else if (second_iterator_type == "U") {
			second_it_min = U;
		}
		else if (second_iterator_type == "V") {
			second_it_min = V;
		}
		else {
			second_it_min = 0;
		}
	}
	void Model::ModelParameters::incrementGlobalIterator()
	{
		incrementer(global_iterator_type, global_step);
		if (second_iterator_type == "T") {
			temperature = second_it_min;
		}
		else if (second_iterator_type == "U") {
			U = second_it_min;
		}
		else if (second_iterator_type == "V") {
			V = second_it_min;
		}
	}
	void Model::ModelParameters::incrementSecondIterator()
	{
		incrementer(second_iterator_type, second_step);
	}
	void Model::ModelParameters::printGlobal() const
	{
		std::cout << global_iterator_type << " = ";
		if (global_iterator_type == "T") {
			std::cout << temperature;
		}
		else if (global_iterator_type == "U") {
			std::cout << U;
		}
		else if (global_iterator_type == "V") {
			std::cout << V;
		}
	}
}