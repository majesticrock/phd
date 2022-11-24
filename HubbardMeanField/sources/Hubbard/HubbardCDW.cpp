#define _USE_MATH_DEFINES

#include "HubbardCDW.hpp"

using Eigen::MatrixXd;

namespace Hubbard {
	HubbardCDW::HubbardCDW(double _temperature, double _U, double _V)
		: BasicHubbardModel(_temperature, _U), V(_V)
	{
		this->delta_cdw = abs(U - V) * 0.5;
		this->delta_sc = abs(U + V) * 0.5;
		if (V > 0) {
			this->delta_sc *= 0.5;
		}
		else {
			this->delta_cdw *= 0.5;
		}
		this->delta_eta = 0.01;

		this->hamilton = MatrixXd::Zero(4, 4);
	}
}
