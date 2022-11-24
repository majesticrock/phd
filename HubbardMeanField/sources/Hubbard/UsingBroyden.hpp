#pragma once
#include "Model.hpp"
#include "Constants.hpp"

namespace Hubbard {
	class UsingBroyden : public Model
	{
	protected:
		double delta_sc, delta_cdw, delta_eta;
		double V;

		double unperturbed_energy(double k_x, double k_y) const override;
		virtual void fillMatrix(double k_x, double k_y) override;
		virtual inline void setParameters(Eigen::VectorXd& F) {
			F(0) *= (this->U - this->V) / (4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);
			F(1) *= (this->U + this->V) / (4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);
			F(2) *= (this->U + this->V) / (4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);
			this->delta_cdw = 0.5 * (F(0) + this->delta_cdw);
			this->delta_sc = 0.5 * (F(1) + this->delta_sc);
			this->delta_eta = 0.5 * (F(2) + this->delta_eta);
		};
	public:
		struct data_set {
			double delta_cdw, delta_sc, delta_eta;
			void print() const;
		};
		UsingBroyden(double _temperature, double _U, double _V);

		data_set compute(const bool print = false);
	};
}
