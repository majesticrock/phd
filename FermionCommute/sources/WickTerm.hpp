#pragma once
#include "Term.hpp"

struct WickOperator {
	std::string type;
	Momentum momentum;
	std::vector<std::string> indizes;

	WickOperator(const std::string& _type, const Momentum& _momentum, const std::vector<std::string>& _indizes = std::vector<std::string>());
	WickOperator(const std::string& _type, const Momentum& _momentum, const std::string& _index);
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
	std::vector<std::pair<Momentum, Momentum>> delta_momenta;
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
};

std::ostream& operator<<(std::ostream& os, const WickOperator& op);
std::ostream& operator<<(std::ostream& os, const std::vector<WickOperator>& ops);
std::ostream& operator<<(std::ostream& os, const WickTerm& term);
std::ostream& operator<<(std::ostream& os, const std::vector<WickTerm>& terms);