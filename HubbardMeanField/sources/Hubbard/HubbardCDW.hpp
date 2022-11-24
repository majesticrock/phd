#pragma once
#include "BasicHubbardModel.hpp"

namespace Hubbard {
	class HubbardCDW : public BasicHubbardModel
	{
	protected:
		double V;

		virtual inline void setParameters(double cdw, double sc, double eta) override {
			this->delta_cdw = ((this->U - this->V) / (4. * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION)) * cdw;
			this->delta_sc = ((this->U + this->V) / (4. * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION)) * sc;
			this->delta_eta = ((this->U + this->V) / (4. * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION)) * eta;
		};

	public:
		HubbardCDW(ModelParameters& _params);
	};
}