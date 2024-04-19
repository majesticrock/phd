#pragma once

namespace Utility::Numerics::Integration {
	template <class UnaryFunction, class RealType, class... Args>
	auto trapezoidal_rule(const UnaryFunction& function, const RealType begin, const RealType end, unsigned long num_steps, Args&&... args) {
		const RealType step = (end - begin) / num_steps;
		auto value = 0.5 * (function(begin, std::forward<Args>(args)...) + function(end, std::forward<Args>(args)...));
		for (int n = 1; n < num_steps; ++n) {
			value += function(n * step, std::forward<Args>(args)...);
		}
		value *= step;
		return value;
	}

	template <class UnaryFunction, class RealType>
	auto trapezoidal_rule(const UnaryFunction& function, const RealType begin, const RealType end, unsigned long num_steps) {
		const RealType step = (end - begin) / num_steps;
		auto value = 0.5 * (function(begin) + function(end));
		for (unsigned long n = 1; n < num_steps; ++n) {
			value += function(n * step);
		}
		value *= step;
		return value;
	}
}