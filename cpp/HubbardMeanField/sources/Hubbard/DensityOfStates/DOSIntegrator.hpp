#pragma once
#include "../GlobalDefinitions.hpp"

namespace Hubbard::DensityOfStates {
	template <class ResultType, class DOS>
	class DOSIntegrator {
	private:
		ResultType result;
		ResultType buffer;
		ResultType second_buffer;

	public:
		// Assumes that the function values are already saved in order to the vector
		// Also assumes that the values for -gamma are the same as for gamma
		template <class VectorType>
		inline const ResultType& integrate_vector_symmetric(const VectorType& function_values) {
			result = static_cast<global_floating_type>(DOS::values[0] * DOS::weights[0]) * function_values[0];

			for (size_t i = 1U; i < DOS::size(); ++i)
			{
				result += static_cast<global_floating_type>(DOS::values[i] * DOS::weights[i]) * function_values[i];
			}
			return result;
		};

		// Assumes that the function values are already saved in order to the vector
		template <class VectorType>
		inline const ResultType& integrate_vector(const VectorType& function_values) {
			result = static_cast<global_floating_type>(DOS::values[0] * DOS::weights[0])
				* (function_values[0] + function_values[DOS::size()]);

			for (size_t i = 1U; i < DOS::size(); ++i)
			{
				result += static_cast<global_floating_type>(DOS::values[i] * DOS::weights[i])
					* (function_values[i] + function_values[i + DOS::size()]);
			}
			return result;
		};

		// This function passes the result of F by reference,
		// i.e. F(gamma, result) and expects F to fill result accordingly.
		// It assumes that F(gamma) is symmetric about 0
		template <class UnaryFunction>
		inline const ResultType& integrate_by_reference_symmetric(const UnaryFunction& F) {
			F(static_cast<global_floating_type>(DOS::abscissa.front()), buffer);
			result = static_cast<global_floating_type>(DOS::values[0] * DOS::weights[0]) * buffer;

			for (size_t i = 1U; i < DOS::size(); ++i)
			{
				F(static_cast<global_floating_type>(DOS::abscissa[i]), buffer);
				result += static_cast<global_floating_type>(DOS::values[i] * DOS::weights[i]) * buffer;
			}
			return result;
		};

		// This function passes the result of F by reference,
		// i.e. F(gamma, result) and expects F to fill result accordingly.
		template <class UnaryFunction>
		inline const ResultType& integrate_by_reference(const UnaryFunction& F) {
			F(static_cast<global_floating_type>(DOS::abscissa.front()), buffer);
			F(static_cast<global_floating_type>(-DOS::abscissa.front()), second_buffer);
			result = static_cast<global_floating_type>(DOS::values[0] * DOS::weights[0]) * (buffer + second_buffer);

			for (size_t i = 1U; i < DOS::size(); ++i)
			{
				F(static_cast<global_floating_type>(DOS::abscissa[i]), buffer);
				F(static_cast<global_floating_type>(-DOS::abscissa[i]), second_buffer);
				result += static_cast<global_floating_type>(DOS::values[i] * DOS::weights[i]) * (buffer + second_buffer);
			}
			return result;
		};

		// This function assumes that F returns its result by value.
		// It assumes that F(gamma) is symmetric about 0
		template <class UnaryFunction>
		inline const ResultType& integrate_by_value_symmetric(const UnaryFunction& F) {
			result = static_cast<global_floating_type>(DOS::values[0] * DOS::weights[0]) *
				F(static_cast<global_floating_type>(DOS::abscissa.front()));

			for (size_t i = 1U; i < DOS::size(); ++i)
			{
				result += static_cast<global_floating_type>(DOS::values[i] * DOS::weights[i]) *
					F(static_cast<global_floating_type>(DOS::abscissa[i]));
			}
			return result;
		};

		// This function assumes that F returns its result by value.
		template <class UnaryFunction>
		inline const ResultType& integrate_by_value(const UnaryFunction& F) {
			result = static_cast<global_floating_type>(DOS::values[0] * DOS::weights[0]) *
				(F(static_cast<global_floating_type>(DOS::abscissa.front())) + F(static_cast<global_floating_type>(-DOS::abscissa.front())));

			for (size_t i = 1U; i < DOS::size(); ++i)
			{
				result += static_cast<global_floating_type>(DOS::values[i] * DOS::weights[i]) *
					(F(static_cast<global_floating_type>(DOS::abscissa[i])) + F(static_cast<global_floating_type>(-DOS::abscissa[i])));
			}
			return result;
		};

		// This functions passes the index i to the corresponding abscissa to F rather than gamma
		// It assumes that F(gamma) is symmetric about 0
		template <class UnaryFunction>
		inline const ResultType& integrate_by_index_symmetric(const UnaryFunction& F) {
			result = static_cast<global_floating_type>(DOS::values[0] * DOS::weights[0]) * F(0U);

			for (size_t i = 1U; i < DOS::size(); ++i)
			{
				result += static_cast<global_floating_type>(DOS::values[i] * DOS::weights[i]) * F(i);
			}
			return result;
		};

		// This functions passes the index i to the corresponding abscissa to F rather than gamma
		// Additionally, it assumes that F returns F(gamma) by value
		template <class UnaryFunction>
		inline const ResultType& integrate_by_index(const UnaryFunction& F) {
			result = static_cast<global_floating_type>(DOS::values[0] * DOS::weights[0]) * F(0U);

			for (size_t i = 1U; i < DOS::size(); ++i)
			{
				result += static_cast<global_floating_type>(DOS::values[i] * DOS::weights[i]) * (F(i) + F(i + DOS::size()));
			}
			return result;
		};

		DOSIntegrator() = default;
		// If it is neccessary for the ResultType to be initaliazed,
		// e.g. give a vector a certain size. ResultType needs to have a copy constructor though
		DOSIntegrator(const ResultType& initialize_result)
			: result(initialize_result), buffer(initialize_result), second_buffer(initialize_result)
		{};
	};
}