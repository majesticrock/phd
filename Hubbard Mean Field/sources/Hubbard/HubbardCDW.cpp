#define _USE_MATH_DEFINES

#include "HubbardCDW.hpp"

using Eigen::MatrixXd;

namespace Hubbard {
	HubbardCDW::HubbardCDW(double _temperature, double _U, double _V)
		: BasicHubbardModel(_temperature, _U), V(_V)
	{
		this->delta_cdw = 0.1;
		this->delta_sc = 0.1;
		this->delta_eta = 0;

		this->hamilton = MatrixXd::Zero(4, 4);
	}
}
