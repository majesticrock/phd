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

	Coefficient() : name(""), momentum(), indizes(), isDaggered(false) {};
	Coefficient(std::string _name, const Momentum& _momentum, const std::vector<std::string>& _indizes, bool _isDaggered)
		: name(_name), momentum(_momentum), indizes(_indizes), isDaggered(_isDaggered) {};
	Coefficient(std::string _name, char _momentum, bool add_Q, const std::vector<std::string>& _indizes, bool _isDaggered)
		: name(_name), momentum(_momentum, add_Q), indizes(_indizes), isDaggered(_isDaggered) { };
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
	Coefficient coefficient;
	std::vector<char> sum_momenta;
	std::vector<std::string> sum_indizes;
	std::vector<Operator> operators;
	// symbolises the Kronecker delta
	std::vector<std::pair<Momentum, Momentum>> delta_momentum;
	std::vector<std::pair<std::string, std::string>> delta_index;

public:
	int multiplicity;
	Term(int _multiplicity, Coefficient _coefficient, const std::vector<char>& _sum_momenta, const std::vector<std::string>& _sum_indizes, const std::vector<Operator>& _operators = std::vector<Operator>())
		: coefficient(_coefficient), sum_momenta(_sum_momenta), sum_indizes(_sum_indizes), operators(_operators), multiplicity(_multiplicity) {};

	Term(int _multiplicity, Coefficient _coefficient, const std::vector<char>& _sum_momenta, const std::vector<Operator>& _operators = std::vector<Operator>())
		: coefficient(_coefficient), sum_momenta(_sum_momenta), operators(_operators), multiplicity(_multiplicity) {};

	Term(int _multiplicity, Coefficient _coefficient, const std::vector<std::string>& _sum_indizes, const std::vector<Operator>& _operators = std::vector<Operator>())
		: coefficient(_coefficient), sum_indizes(_sum_indizes), operators(_operators), multiplicity(_multiplicity) {};

	Term(int _multiplicity, Coefficient _coefficient, const std::vector<Operator>& _operators = std::vector<Operator>())
		: coefficient(_coefficient), operators(_operators), multiplicity(_multiplicity) {};

	Term() : coefficient(), operators(), multiplicity(0) {};

	inline bool isIdentity() const {
		return this->operators.empty();
	}

	void print() const;
	inline void flipSign() {
		this->multiplicity *= -1;
	}

	void setDeltas();
	void sort();
	// Checks for equality of everything except of multiplicity
	inline bool isEqual(const Term& other) const {
		if (this->coefficient != other.coefficient) return false;
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

		if (this->delta_index.size() != other.delta_index.size()) return false;
		if (this->delta_momentum.size() != other.delta_momentum.size()) return false;
		for (size_t i = 0; i < delta_index.size(); i++)
		{
			if (this->delta_index[i] != other.delta_index[i]) return false;
		}
		for (size_t i = 0; i < delta_momentum.size(); i++)
		{
			if (this->delta_momentum[i] != other.delta_momentum[i]) return false;
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
	friend void commutator(std::vector<Term>& reciever, const std::vector<Term>& left, const std::vector<Term>& right);
	friend std::ostream& operator<<(std::ostream& os, const Term& term);
};

inline bool operator==(const Term& lhs, const Term& rhs) {
	return lhs.isEqual(rhs);
};
inline bool operator!=(const Term& lhs, const Term& rhs) {
	return !(lhs.isEqual(rhs));
};

std::ostream& operator<<(std::ostream& os, const Coefficient& coeff);
std::ostream& operator<<(std::ostream& os, const std::vector<Term>& terms);

void cleanUp(std::vector<Term>& terms);