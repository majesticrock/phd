#pragma once
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>

namespace Utility {
	template <typename T>
	std::string to_string_with_precision(const T a_value, const int n)
	{
		std::ostringstream out;
		out.precision(n);
		out << std::fixed << a_value;
		return out.str();
	}

	void saveData(std::vector<double>& data, std::string filename);
	void saveData(std::vector<std::vector<double>>& data, std::string filename);
	void saveData(std::vector<std::vector<double>>& data, std::string filename, std::vector<std::string> comments);
	void saveData(std::vector<double>& data, const int linebreak, std::string filename, std::vector<std::string> comments);
}