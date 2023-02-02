#pragma once
#include "Model.hpp"
#include "Constants.hpp"

namespace Hubbard {
	class UsingBroyden : public Model
	{
	private:
		Eigen::Vector4d new_sc;
		Eigen::Vector4d new_eta;
		Eigen::Vector4d old_sc;
		Eigen::Vector4d old_eta;
	protected:
		double V;

		virtual void fillHamiltonian(double k_x, double k_y) override;
		virtual inline void setParameters(Eigen::VectorXd& F) {
			F(0) *= (this->U - this->V) / (8 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);
			F(1) *= (this->U) / (8 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION); // + this->V
			F(2) *= (this->U) / (8 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION); // + this->V
			this->delta_cdw = 0.5 * (F(0) + this->delta_cdw);
			this->delta_sc = 0.5 * (F(1) + this->delta_sc);
			this->delta_eta = 0.5 * (F(2) + this->delta_eta);
		};

		virtual void compute_quartics() override;
		virtual void fill_M_N() override;
	public:
		UsingBroyden(ModelParameters& _params);

		data_set computePhases(const bool print = false) override;
	};
}