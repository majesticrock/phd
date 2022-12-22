#pragma once
#include "Model.hpp"
#include "Constants.hpp"

namespace Hubbard {
	class BasicHubbardModel : public Model
	{
	protected:
		double unperturbed_energy(double k_x, double k_y) const override;
		virtual void fillMatrix(double k_x, double k_y) override;
		virtual inline void setParameters(double cdw, double sc, double eta) {
			this->delta_cdw = cdw * this->U / (4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);
			this->delta_sc = sc * this->U / (4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);
			this->delta_eta = eta * this->U / (4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);
		};
	public:
		BasicHubbardModel(ModelParameters& _params);

		data_set compute(const bool print = false);
	};
}