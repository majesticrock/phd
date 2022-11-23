#include "OutputWriter.hpp"
#include <fstream>
#include <iostream>
#include <ctime>

namespace Utility {
	inline std::tm localtime_xp(std::time_t timer)
	{
		std::tm bt{};
#if defined(__unix__)
		localtime_r(&timer, &bt);
#elif defined(_MSC_VER)
		localtime_s(&bt, &timer);
#else
		static std::mutex mtx;
		std::lock_guard<std::mutex> lock(mtx);
		bt = *std::localtime(&timer);
#endif
		return bt;
	}

	// default = "YYYY-MM-DD HH:MM:SS"
	inline std::string time_stamp(const std::string& fmt = "%F %T")
	{
		auto bt = localtime_xp(std::time(0));
		char buf[64];
		return { buf, std::strftime(buf, sizeof(buf), fmt.c_str(), &bt) };
	}

	void saveData(std::vector<double>& data, std::string filename)
	{
		std::ofstream out(filename);
		if (out.is_open()) {
			out << "# " << time_stamp() << "\n\n";

			out << std::scientific << std::setprecision(10);
			for (int i = 0; i < data.size(); i++)
			{
				out << data[i];
				out << " ";
			}
		}
		else {
			std::cerr << "Could not open output filestream for file: " << filename << std::endl;
		}
	}
	void saveData(std::vector<std::vector<double>>& data, std::string filename)
	{
		std::ofstream out(filename);
		if (out.is_open()) {
			out << "# " << time_stamp() << "\n\n";

			out << std::scientific << std::setprecision(10);
			for (int i = 0; i < data.size(); i++)
			{
				for (int s = 0; s < data[0].size(); s++)
				{
					out << data[i][s];
					if (s < data[0].size() - 1) {
						out << " ";
					}
				}
				out << "\n";
			}
		}
		else {
			std::cerr << "Could not open output filestream for file: " << filename << std::endl;
		}
	}
	void saveData(std::vector<std::vector<double>>& data, std::string filename, std::vector<std::string> comments)
	{
		std::ofstream out(filename);
		if (out.is_open()) {
			out << "# " << time_stamp() << "\n#\n";

			for (int i = 0; i < comments.size(); i++)
			{
				out << "# " << comments[i] << "\n";
			}

			out << "\n" << std::scientific << std::setprecision(10);
			for (int i = 0; i < data.size(); i++)
			{
				for (int s = 0; s < data[0].size(); s++)
				{
					out << data[i][s];
					if (s < data[0].size() - 1) {
						out << " ";
					}
				}
				out << "\n";
			}
		}
		else {
			std::cerr << "Could not open output filestream for file: " << filename << std::endl;
		}
	}
	void saveData(std::vector<double>& data, const int linebreak, std::string filename, std::vector<std::string> comments) {
		std::ofstream out(filename);
		if (out.is_open()) {
			out << "# " << time_stamp() << "\n#\n";

			for (int i = 0; i < comments.size(); i++)
			{
				out << "# " << comments[i] << "\n";
			}

			out << std::scientific << std::setprecision(10);
			for (int i = 0; i < data.size(); i++)
			{
				if(i % linebreak == 0){
					out << "\n";
				}
				out << data[i];
				out << " ";
			}
		}
		else {
			std::cerr << "Could not open output filestream for file: " << filename << std::endl;
		}
	}
}