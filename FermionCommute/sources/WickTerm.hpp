#pragma once
#include "Term.hpp"

typedef std::pair<Momentum, Momentum> pair_of_momenta;
struct WickOperator {
	std::string type;
	bool isDaggered;
	Momentum momentum;
	std::vector<std::string> indizes;

	WickOperator(const std::string& _type, const bool _isDaggered, const Momentum& _momentum, const std::vector<std::string>& _indizes = std::vector<std::string>());
	WickOperator(const std::string& _type, const bool _isDaggered, const Momentum& _momentum, const std::string& _index);
	WickOperator();
};

class Term;

struct WickTerm
{
	int multiplicity;
	std::vector<Coefficient> coefficients;
	std::vector<char> sum_momenta;
	std::vector<std::string> sum_indizes;
	std::vector<WickOperator> operators;

	// symbolises the Kronecker delta
	std::vector<pair_of_momenta> delta_momenta;
	std::vector<std::pair<std::string, std::string>> delta_indizes;

	std::vector<Operator> temporary_operators;

	explicit WickTerm(const Term* base);

	inline bool isIdentity() const {
		return this->operators.empty();
	}
	inline bool handled() const {
		if (this->temporary_operators.empty()) return true;
		return !(this->operators.empty());
	}
	bool swapToWickOperators(std::vector<WickTerm>& reciever);
	bool setDeltas();
	void computeSums();
	void discardZeroMomenta();
	void renameSums();
	void sort();
};

inline bool operator==(const WickOperator& lhs, const WickOperator& rhs) {
	if (lhs.type != rhs.type) return false;
	if (lhs.isDaggered != rhs.isDaggered) return false;
	if (lhs.momentum != rhs.momentum) return false;
	if (lhs.indizes.size() != rhs.indizes.size()) return false;
	for (size_t i = 0; i < lhs.indizes.size(); i++)
	{
		if (lhs.indizes[i] != rhs.indizes[i]) return false;
	}
	return true;
};
inline bool operator!=(const WickOperator& lhs, const WickOperator& rhs) {
	return !(lhs == rhs);
};
inline bool operator==(const WickTerm& lhs, const WickTerm& rhs) {
	if (lhs.coefficients.size() != rhs.coefficients.size()) return false;
	for (size_t i = 0; i < lhs.coefficients.size(); i++)
	{
		if (lhs.coefficients[i] != rhs.coefficients[i]) return false;
	}
	if (lhs.sum_indizes.size() != rhs.sum_indizes.size()) return false;
	if (lhs.sum_momenta.size() != rhs.sum_momenta.size()) return false;
	for (size_t i = 0; i < lhs.sum_indizes.size(); i++)
	{
		if (lhs.sum_indizes[i] != rhs.sum_indizes[i]) return false;
	}
	for (size_t i = 0; i < lhs.sum_momenta.size(); i++)
	{
		if (lhs.sum_momenta[i] != rhs.sum_momenta[i]) return false;
	}

	if (lhs.delta_indizes.size() != rhs.delta_indizes.size()) return false;
	if (lhs.delta_momenta.size() != rhs.delta_momenta.size()) return false;
	for (size_t i = 0; i < lhs.delta_indizes.size(); i++)
	{
		if (lhs.delta_indizes[i] != rhs.delta_indizes[i]) return false;
	}
	for (size_t i = 0; i < lhs.delta_momenta.size(); i++)
	{
		if (lhs.delta_momenta[i] != rhs.delta_momenta[i]) return false;
	}

	if (lhs.operators.size() != rhs.operators.size()) return false;

	// The Wick "operators" are actually just numbers
	// therefore I might be interested to implement permutations as well...
	for (size_t i = 0; i < lhs.operators.size(); i++)
	{
		if (lhs.operators[i] != rhs.operators[i]) return false;
	}
	return true;
};
inline bool operator!=(const WickTerm& lhs, const WickTerm& rhs) {
	return !(lhs == rhs);
};

void cleanWicks(std::vector<WickTerm>& terms);

std::ostream& operator<<(std::ostream& os, const WickOperator& op);
std::ostream& operator<<(std::ostream& os, const std::vector<WickOperator>& ops);
std::ostream& operator<<(std::ostream& os, const WickTerm& term);
std::ostream& operator<<(std::ostream& os, const std::vector<WickTerm>& terms);