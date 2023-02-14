#pragma once
#include <string>

#include "Operator.hpp"

struct Coefficient {
	std::string name;
	Momentum momentum;
	// Contains all indizes, standard: first index = spin, all others arbitrary, e.g. orbitals, bands etc
	std::vector<std::string> indizes;
	bool isDaggered;
	// if Coeff(k) = Coeff(-k)
	bool translationalInvariance = true;

	Coefficient();
	explicit Coefficient(std::string _name);
	Coefficient(std::string _name, const Momentum& _momentum, const std::vector<std::string>& _indizes, bool _isDaggered = false);
	Coefficient(std::string _name, char _momentum, bool add_Q, const std::vector<std::string>& _indizes, bool _isDaggered = false);
	Coefficient(std::string _name, const Momentum& _momentum, bool _isDaggered = false);
	Coefficient(std::string _name, char _momentum, bool add_Q = false, bool _isDaggered = false);
};

inline bool operator==(const Coefficient& lhs, const Coefficient& rhs) {
	if (lhs.name != rhs.name) return false;
	if (lhs.momentum != rhs.momentum) return false;
	if (lhs.isDaggered != rhs.isDaggered) return false;
	for (size_t i = 0; i < lhs.indizes.size(); i++)
	{
		if (lhs.indizes[i] != rhs.indizes[i]) return false;
	}
	return true;
}
inline bool operator!=(const Coefficient& lhs, const Coefficient& rhs) {
	return !(lhs == rhs);
}