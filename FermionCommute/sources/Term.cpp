#define append_vector(a, b) a.insert(a.end(), b.begin(), b.end())
#include "Term.hpp"
#include <cmath>

void Term::print() const {
	std::cout << *this << std::endl;
}

void Term::setDeltas()
{
	// Set all deltas up to the same notation
	for (auto& delta : delta_momentum) {
		if (delta.first.add_Q) {
			delta.first.add_Q = false;
			delta.second.add_Q = (delta.second.add_Q != true);
		}
		if (delta.second.momentum_list.size() == 1 && delta.first.momentum_list.size() > 1) {
			std::swap(delta.first, delta.second);
		}
		if (delta.first.momentum_list.size() == 1 && delta.first.momentum_list[0].first < 0) {
			delta.first.flipMomentum();
			delta.second.flipMomentum();
		}
		else if (delta.first.momentum_list.size() > 1 && delta.second.momentum_list.size() > 1) {
			bool foundCandidate = false;
			int index = 0;
			delta.second -= delta.first;
			delta.first.momentum_list.clear();

			for (auto m : sum_momenta)
			{
				index = delta.second.isUsed(m);
				if (index >= 0) {
					foundCandidate = true;
					if (abs(delta.second.momentum_list[index].first) == 1) {
						break;
					}
				}
				if (!foundCandidate) {
					index = 0;
				}
				if (delta.second.momentum_list[index].first > 0) {
					delta.second.flipMomentum();
				}
				delta.first.momentum_list.push_back(delta.second.momentum_list[index]);
				delta.first.flipMomentum();
				delta.second.momentum_list.erase(delta.second.momentum_list.begin() + index);
			}
		}

		for (auto& op : operators) {
			op.momentum.replaceOccurances(delta.first.momentum_list[0].second, delta.second);
		}
	}
	for (auto& delta : delta_index) {
		for (auto& op : operators) {
			for (auto it = op.indizes.begin(); it != op.indizes.end(); ++it)
			{
				if (delta.first == UP || delta.first == DOWN) {
					if (*it == delta.second) {
						*it = delta.first;
					}
				}
				else {
					if (*it == delta.first) {
						*it = delta.second;
					}
				}
			}
		}
	}

	// Remove delta^2
	for (size_t i = 0; i < delta_momentum.size(); i++)
	{
	setDeltas_outerLoop_momentum:

		for (size_t j = i + 1; j < delta_momentum.size(); j++)
		{
			if (delta_momentum[i] == delta_momentum[j]) {
				delta_momentum.erase(delta_momentum.begin() + j);
				goto setDeltas_outerLoop_momentum;
			}
		}
	}
	for (size_t i = 0; i < delta_index.size(); i++)
	{
	setDeltas_outerLoop_index:

		for (size_t j = i + 1; j < delta_index.size(); j++)
		{
			if (delta_index[i] == delta_index[j]) {
				delta_index.erase(delta_index.begin() + j);
				goto setDeltas_outerLoop_index;
			}
		}
	}
}

void Term::sort()
{
	for (int i = 0; i < operators.size(); i++)
	{
		for (int j = i + 1; j < operators.size(); j++)
		{
			if (operators[i].isDaggered == operators[j].isDaggered) {
				if (operators[i].isDaggered) {
					// c^+ c^+
					if (operators[j].indizes[0] == UP && operators[i].indizes[0] != UP) {
						std::swap(operators[i], operators[j]);
						if (abs(i - j) % 2 != 0) {
							flipSign();
						}
					}
				}
				else {
					// c c
					if (operators[j].indizes[0] == DOWN && operators[i].indizes[0] != DOWN) {
						std::swap(operators[i], operators[j]);
						if (abs(i - j) % 2 != 0) {
							flipSign();
						}
					}
				}
			}
		}
	}
}

void normalOrder(std::vector<Term>& terms) {
	for (size_t t = 0; t < terms.size();) {
	normalOder_outerLoop:
		size_t n = terms[t].operators.size();
		size_t new_n;
		while (n > 1) {
			new_n = 0;
			for (size_t i = 1; i < terms[t].operators.size(); i++)
			{
				if (!(terms[t].operators[i - 1].isDaggered) && (terms[t].operators[i].isDaggered)) {
					bool other_deltas = false;
					new_n = i;
					// Swap cc^+
					terms[t].flipSign();
					std::swap(terms[t].operators[i - 1], terms[t].operators[i]);

					// Add a new term where cc^+ is replaced by the appropriate delta
					Term new_term(terms[t]);
					new_term.flipSign();
					if (new_term.operators[i - 1].indizes.size() != new_term.operators[i].indizes.size()) {
						std::cerr << "Operators do not have the same index count." << std::endl;
						throw;
					}

					if ((new_term.operators[i - 1].indizes[0] == UP || new_term.operators[i - 1].indizes[0] == DOWN)
						&& (new_term.operators[i].indizes[0] == UP || new_term.operators[i].indizes[0] == DOWN)) {
						if (new_term.operators[i - 1].indizes[0] != new_term.operators[i].indizes[0]) {
							// Case: one up the other down, then free anticommutation
							continue;
						}
					}
					else {
						new_term.delta_index.push_back(
							std::make_pair(new_term.operators[i - 1].indizes[0], new_term.operators[i].indizes[0]));
					}
					for (size_t c = 1; c < new_term.operators[i - 1].indizes.size(); c++)
					{
						// if the indizes are not the same we emplace a delta
						// otherwise no action is required
						if (new_term.operators[i - 1].indizes[c] != new_term.operators[i].indizes[c]) {
							other_deltas = true;
							new_term.delta_index.push_back(
								std::make_pair(new_term.operators[i - 1].indizes[c], new_term.operators[i].indizes[c]));
						}
					}
					if (new_term.operators[i - 1].momentum != new_term.operators[i].momentum) {
						other_deltas = true;
						new_term.delta_momentum.push_back(
							std::make_pair(new_term.operators[i - 1].momentum, new_term.operators[i].momentum)
						);
					}

					new_term.operators.erase(new_term.operators.begin() + i - 1, new_term.operators.begin() + i + 1);
					if (other_deltas) terms.push_back(new_term);
				}
				else if (terms[t].operators[i - 1] == terms[t].operators[i]) {
					// two identical fermion operators = 0
					terms.erase(terms.begin() + t);
					goto normalOder_outerLoop;
				}
			}
			n = new_n;
		}
		++t;
	}
}

void commutator(std::vector<Term>& reciever, const Term& left, const Term& right)
{
	reciever.resize(2);
	reciever[0] = left;
	append_vector(reciever[0].operators, right.operators);
	reciever[1] = right;
	append_vector(reciever[1].operators, left.operators);
	reciever[1].flipSign();
	normalOrder(reciever);
}

void commutator(std::vector<Term>& reciever, const std::vector<Term>& left, const std::vector<Term>& right)
{
	const size_t HALF_SIZE = left.size() * right.size();
	reciever.resize(2 * HALF_SIZE);
	for (size_t i = 0; i < left.size(); i++)
	{
		for (size_t j = 0; j < right.size(); j++)
		{
			reciever[i * right.size() + j] = left[i];
			append_vector(reciever[i * right.size() + j].operators, right[j].operators);

			reciever[HALF_SIZE + i * right.size() + j] = right[j];
			append_vector(reciever[HALF_SIZE + i * right.size() + j].operators, left[i].operators);
			reciever[HALF_SIZE + i * right.size() + j].flipSign();
		}
	}
	normalOrder(reciever);
}

std::ostream& operator<<(std::ostream& os, const Term& term)
{
	if (term.multiplicity > 0) {
		os << "+";
	}
	os << term.multiplicity << " \\cdot ";
	if (!term.sum_indizes.empty()) {
		os << "\\sum_{ ";
		for (const auto& index : term.sum_indizes) {
			os << index << " ";
		}
		os << "}";
	}
	if (!term.sum_momenta.empty()) {
		os << "\\sum_{ ";
		for (const auto& momentum : term.sum_momenta) {
			os << momentum << " ";
		}
		os << "}";
	}
	os << term.coefficient << " ";
	for (const auto& delta : term.delta_momentum) {
		os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
	}
	for (const auto& delta : term.delta_index) {
		os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
	}

	if (term.isIdentity()) {
		os << " \\mathbb{1} ";
		return os;
	}
	for (const auto& op : term.operators) {
		os << op << " ";
	}
	return os;
}

std::ostream& operator<<(std::ostream& os, const Coefficient& coeff)
{
	os << coeff.name;
	if (!coeff.indizes.empty()) {
		os << "_{ ";
		for (const auto& index : coeff.indizes) {
			os << index << " ";
		}
		os << "}";
	}
	if (coeff.isDaggered) {
		os << "^*";
	}
	if (!coeff.momentum.momentum_list.empty()) {
		os << " ( " << coeff.momentum << " )";
	}
	return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<Term>& terms)
{
	for (std::vector<Term>::const_iterator it = terms.begin(); it != terms.end(); ++it)
	{
		os << "\t&" << *it;
		if (it != terms.end() - 1) {
			os << " \\\\";
		}
		os << "\n";
	}
	return os;
}

void cleanUp(std::vector<Term>& terms)
{
	for (auto& term : terms) {
		term.setDeltas();
		term.sort();
	}

	for (size_t i = 0; i < terms.size(); i++)
	{
	cleanUp_outerLoop:

		for (size_t j = i + 1; j < terms.size(); j++)
		{
			if (terms[i] == terms[j]) {
				terms[i].multiplicity += terms[j].multiplicity;
				terms.erase(terms.begin() + j);
				goto cleanUp_outerLoop;
			}
		}
	}
}