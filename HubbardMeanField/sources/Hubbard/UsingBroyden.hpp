#pragma once
#include "Model.hpp"
#include "Constants.hpp"

namespace Hubbard {
	class UsingBroyden : public Model
	{
	protected:
		double V;

		double unperturbed_energy(double k_x, double k_y) const override;
		virtual void fillHamiltonian(double k_x, double k_y) override;
		virtual inline void setParameters(Eigen::VectorXd& F) {
			F(0) *= (this->U - this->V) / (8 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);
			F(1) *= (this->U + 0.5 * this->V) / (8 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);
			F(2) *= (this->U + 0.5 * this->V) / (8 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);
			this->delta_cdw = 0.5 * (F(0) + this->delta_cdw);
			this->delta_sc = 0.5 * (F(1) + this->delta_sc);
			this->delta_eta = 0.5 * (F(2) + this->delta_eta);
		};
	public:
		UsingBroyden(ModelParameters& _params);

		data_set computePhases(const bool print = false) override;
	};
}