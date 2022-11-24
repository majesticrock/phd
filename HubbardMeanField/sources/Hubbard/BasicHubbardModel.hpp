#pragma once
#include "Model.hpp"
#include "Constants.hpp"

namespace Hubbard {
	class BasicHubbardModel : public Model
	{
	protected:
		double delta_sc, delta_cdw, delta_eta;

		double unperturbed_energy(double k_x, double k_y) const override;
		virtual void fillMatrix(double k_x, double k_y) override;
		virtual inline void setParameters(double cdw, double sc, double eta) {
			this->delta_cdw = cdw * this->U / (4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);
			this->delta_sc = sc * this->U / (4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);
			this->delta_eta = eta * this->U / (4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);
		};
	public:
		struct data_set {
			double delta_cdw, delta_sc, delta_eta;
			void print() const;
		};
		BasicHubbardModel(ModelParameters& _params);

		data_set compute(const bool print = false);
	};
}