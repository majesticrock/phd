#define _USE_MATH_DEFINES

#include "HubbardCDW.hpp"

using Eigen::MatrixXd;

namespace Hubbard {
	HubbardCDW::HubbardCDW(double _temperature, double _U, double _V)
		: BasicHubbardModel(_temperature, _U), V(_V)
	{
		this->delta_sc  = (abs(abs(this->U) - this->temperature) + 0.1) * 0.05;
		this->delta_cdw = (abs(abs(this->U) - this->temperature) + 0.1) * 0.05;
		this->delta_eta = (abs(abs(this->U) - this->temperature) + 0.1) * 0.01;

		this->hamilton = MatrixXd::Zero(4, 4);
	}
}
