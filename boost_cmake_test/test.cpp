#include <boost/archive/binary_oarchive.hpp>
#include <iostream>
#include <fstream>

int main(int argc, char** argv) {
    std::ofstream ofs("build/test.bin", std::ios::binary);
	boost::archive::binary_oarchive oa(ofs);
	oa << "Just some test string";
	ofs.close();

	std::cout << "Succesfully executed the program!" << std::endl;

    return 0;
}