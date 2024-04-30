#pragma once
#include "_internal_functions.hpp"
#include "../Resolvent.hpp"

namespace Utility::Numerics::iEoM {
	template <class NumberType>
	class IEoMResolventInterface {
	private:
		using Matrix = Eigen::Matrix<NumberType, Eigen::Dynamic, Eigen::Dynamic>;
		using Vector = Eigen::Vector<NumberType, Eigen::Dynamic>;
		using RealType = UnderlyingFloatingPoint_t<NumberType>;
	protected:
		ieom_internal<NumberType> _internal;

		virtual void fill_N() = 0;
		virtual void fill_M() = 0;
		virtual void createStartingStates() = 0;

		inline void fillMatrices() {
			this->fill_N();
			this->fill_M();
		}
	public:
		virtual bool dynamic_matrix_is_negative() = 0;
		virtual std::vector<ResolventDataWrapper<RealType>> computeCollectiveModes(unsigned int LANCZOS_ITERATION_NUMBER) = 0;

		IEoMResolventInterface(RealType const& sqrt_precision)
			: _internal(sqrt_precision) {};
	};
}