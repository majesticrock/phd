#pragma once
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include <limits>
#include "../UnderlyingFloatingPoint.hpp"

namespace Utility::Selfconsistency {
    template<class RealType>
    constexpr RealType PRECISION = 1e2 * std::numeric_limits<RealType>::epsilon();
    template<> constexpr float PRECISION<float> = 1e1 * std::numeric_limits<float>::epsilon();

    struct DebugPolicy {
        bool convergenceWarning{true};
        bool printSteps{false};
    };
    constexpr DebugPolicy WarnNoConvergence{ false, true };
	constexpr DebugPolicy NoWarning{ false, false };
	constexpr DebugPolicy PrintEverything{ true, true };

	template <typename DataType, class Model, class SelfconsistencyAttributes, const DebugPolicy& debugPolicy=WarnNoConvergence>
	class IterativeSolver {
	protected:
		Model* _model{};
		SelfconsistencyAttributes* _attr{};
		const size_t NUMBER_OF_PARAMETERS;

		using ParameterVector = Eigen::Vector<DataType, Eigen::Dynamic>;
		using RealType = UnderlyingFloatingPoint_t<DataType>;

		inline bool hasSignFlippingBehaviour(const ParameterVector& x0) {
			for (size_t j = 0U; j < this->NUMBER_OF_PARAMETERS; ++j)
			{
				if (abs(x0[j]) > 1e1 * PRECISION<DataType>) {
					if (abs((x0[j] + (*this->_attr)[j]) / x0[j]) < PRECISION<DataType>) {
						return true;
					}
				}
			}
			return false;
		};
		bool procedureIterative(const size_t MAX_STEPS)
		{
			UnderlyingFloatingPoint_t<DataType> error{ 100 };

			ParameterVector f0{ ParameterVector::Zero(this->NUMBER_OF_PARAMETERS) };
			std::copy(this->_attr->begin(), this->_attr->end(), f0.begin());
			ParameterVector x0{ f0 };

			if (debugPolicy.printAll) {
				std::cout << "-1:\t" << std::scientific << std::setprecision(4) << "\n" << x0.transpose() << std::endl;
			}

			size_t iterNum = 0U;
			while (iterNum < MAX_STEPS && error > PRECISION<DataType>) {
				this->_model->iterationStep(x0, f0);
				if (hasSignFlippingBehaviour(x0)) {
					if constexpr (debugPolicy.convergenceWarning) {
						std::cerr << "Sign flipper for " << this->_model->info() << std::endl;
					}
					this->_attr->setZero();
					return false;
				}

				error = f0.norm();
				std::copy(this->_attr->begin(), this->_attr->end(), x0.begin());

				if constexpr (debugPolicy.printAll) {
					std::cout << iterNum << ":\t" << std::scientific << std::setprecision(4) << "\n" << x0.transpose() << std::endl;
					std::cout << "Error:  " << error << "\n";
				}
				++iterNum;
			}
			if (iterNum >= MAX_STEPS) {
				return false;
			}
			return true;
		};
	public:
		virtual const SelfconsistencyAttributes& computePhases()
		{
			constexpr size_t MAX_STEPS = 1500;
			this->_attr->converged = true;
			if (!procedureIterative(MAX_STEPS)) {
				if constexpr (debugPolicy.convergenceWarning) {
					std::cerr << "No convergence for " << this->_model->info() << std::endl;
				}
				this->_attr->setZero();
			}

			return *this->_attr;
		};

		virtual ~IterativeSolver() = default;
		IterativeSolver() = delete;
		IterativeSolver(Model* model_ptr, SelfconsistencyAttributes* attribute_ptr)
			: _model(model_ptr), _attr(attribute_ptr), NUMBER_OF_PARAMETERS(_attr->size()) {};
	};
}