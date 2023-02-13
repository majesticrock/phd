#pragma once
#include <vector>
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

	Coefficient() : name(""), momentum(), indizes(), isDaggered(false) {};

	explicit Coefficient(std::string _name) : name(_name), momentum(), indizes(), isDaggered(false) {};

	Coefficient(std::string _name, const Momentum& _momentum, const std::vector<std::string>& _indizes, bool _isDaggered = false)
		: name(_name), momentum(_momentum), indizes(_indizes), isDaggered(_isDaggered) {};

	Coefficient(std::string _name, char _momentum, bool add_Q, const std::vector<std::string>& _indizes, bool _isDaggered = false)
		: name(_name), momentum(_momentum, add_Q), indizes(_indizes), isDaggered(_isDaggered) { };

	Coefficient(std::string _name, const Momentum& _momentum, bool _isDaggered = false)
		: name(_name), momentum(_momentum), indizes(), isDaggered(_isDaggered) {};

	Coefficient(std::string _name, char _momentum, bool add_Q = false, bool _isDaggered = false)
		: name(_name), momentum(_momentum, 1, add_Q), indizes(), isDaggered(_isDaggered) { };
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

class Term {
private:
	std::vector<Coefficient> coefficients;
	std::vector<char> sum_momenta;
	std::vector<std::string> sum_indizes;
	std::vector<Operator> operators;
	// symbolises the Kronecker delta
	std::vector<std::pair<Momentum, Momentum>> delta_momenta;
	std::vector<std::pair<std::string, std::string>> delta_indizes;

public:
	int multiplicity;
	Term(int _multiplicity, Coefficient _coefficient, const std::vector<char>& _sum_momenta, const std::vector<std::string>& _sum_indizes, const std::vector<Operator>& _operators = std::vector<Operator>())
		: coefficients(1, _coefficient), sum_momenta(_sum_momenta), sum_indizes(_sum_indizes), operators(_operators), multiplicity(_multiplicity) {};

	Term(int _multiplicity, Coefficient _coefficient, const std::vector<char>& _sum_momenta, const std::vector<Operator>& _operators = std::vector<Operator>())
		: coefficients(1, _coefficient), sum_momenta(_sum_momenta), operators(_operators), multiplicity(_multiplicity) {};

	Term(int _multiplicity, Coefficient _coefficient, const std::vector<std::string>& _sum_indizes, const std::vector<Operator>& _operators = std::vector<Operator>())
		: coefficients(1, _coefficient), sum_indizes(_sum_indizes), operators(_operators), multiplicity(_multiplicity) {};

	Term(int _multiplicity, Coefficient _coefficient, const std::vector<Operator>& _operators = std::vector<Operator>())
		: coefficients(1, _coefficient), operators(_operators), multiplicity(_multiplicity) {};

	Term(int _multiplicity, const std::vector<char>& _sum_momenta, const std::vector<std::string>& _sum_indizes, const std::vector<Operator>& _operators = std::vector<Operator>())
		: coefficients(), sum_momenta(_sum_momenta), sum_indizes(_sum_indizes), operators(_operators), multiplicity(_multiplicity) {};

	Term(int _multiplicity, const std::vector<char>& _sum_momenta, const std::vector<Operator>& _operators = std::vector<Operator>())
		: coefficients(), sum_momenta(_sum_momenta), operators(_operators), multiplicity(_multiplicity) {};

	Term(int _multiplicity, const std::vector<std::string>& _sum_indizes, const std::vector<Operator>& _operators = std::vector<Operator>())
		: coefficients(), sum_indizes(_sum_indizes), operators(_operators), multiplicity(_multiplicity) {};

	explicit Term(int _multiplicity, const std::vector<Operator>& _operators = std::vector<Operator>())
		: coefficients(), operators(_operators), multiplicity(_multiplicity) {};

	Term() : coefficients(), operators(), multiplicity(0) {};

	inline bool isIdentity() const {
		return this->operators.empty();
	}

	void print() const;
	inline void flipSign() {
		this->multiplicity *= -1;
	}

	void setDeltas();
	void computeSums();
	void discardZeroMomenta();
	void sort();
	// Checks for equality of everything except of multiplicity
	inline bool isEqual(const Term& other) const {
		if (this->coefficients.size() != other.coefficients.size()) return false;
		for (size_t i = 0; i < coefficients.size(); i++)
		{
			if (this->coefficients[i] != other.coefficients[i]) return false;
		}
		if (this->sum_indizes.size() != other.sum_indizes.size()) return false;
		if (this->sum_momenta.size() != other.sum_momenta.size()) return false;
		for (size_t i = 0; i < sum_indizes.size(); i++)
		{
			if (this->sum_indizes[i] != other.sum_indizes[i]) return false;
		}
		for (size_t i = 0; i < sum_momenta.size(); i++)
		{
			if (this->sum_momenta[i] != other.sum_momenta[i]) return false;
		}

		if (this->delta_indizes.size() != other.delta_indizes.size()) return false;
		if (this->delta_momenta.size() != other.delta_momenta.size()) return false;
		for (size_t i = 0; i < delta_indizes.size(); i++)
		{
			if (this->delta_indizes[i] != other.delta_indizes[i]) return false;
		}
		for (size_t i = 0; i < delta_momenta.size(); i++)
		{
			if (this->delta_momenta[i] != other.delta_momenta[i]) return false;
		}

		if (this->operators.size() != other.operators.size()) return false;
		for (size_t i = 0; i < this->operators.size(); i++)
		{
			if (this->operators[i] != other.operators[i]) return false;
		}
		return true;
	};

	std::string toStringWithoutPrefactor() const;

	friend void normalOrder(std::vector<Term>& terms);
	friend void commutator(std::vector<Term>& reciever, const Term& left, const Term& right);
	friend std::ostream& operator<<(std::ostream& os, const Term& term);
};
void commutator(std::vector<Term>& reciever, const std::vector<Term>& left, const std::vector<Term>& right);
inline void commutator(std::vector<Term>& reciever, const Term& left, const std::vector<Term>& right) {
	const std::vector<Term> buffer = { left };
	commutator(reciever, buffer, right);
};
inline void commutator(std::vector<Term>& reciever, const std::vector<Term>& left, const Term& right) {
	const std::vector<Term> buffer = { right };
	commutator(reciever, left, buffer);
};

inline bool operator==(const Term& lhs, const Term& rhs) {
	return lhs.isEqual(rhs);
};
inline bool operator!=(const Term& lhs, const Term& rhs) {
	return !(lhs.isEqual(rhs));
};

std::ostream& operator<<(std::ostream& os, const Coefficient& coeff);
std::ostream& operator<<(std::ostream& os, const std::vector<Coefficient>& coeffs);
std::ostream& operator<<(std::ostream& os, const std::vector<Term>& terms);

void cleanUp(std::vector<Term>& terms);