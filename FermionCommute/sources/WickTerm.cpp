#include "WickTerm.hpp"
#include <map>

#define LEFT temporary_operators[i]
#define RIGHT temporary_operators[i + 1]
#define L_SPIN temporary_operators[i].indizes[0]
#define R_SPIN temporary_operators[i + 1].indizes[0]

WickTerm::WickTerm(const Term* base)
	: multiplicity(base->multiplicity), coefficients(base->coefficients), sum_momenta(base->sum_momenta), sum_indizes(base->sum_indizes),
	operators(), delta_momenta(base->delta_momenta), delta_indizes(base->delta_indizes), temporary_operators()
{
}

bool WickTerm::swapToWickOperators(std::vector<WickTerm>& reciever)
{
	this->operators.reserve(temporary_operators.size() / 2);
	WickTerm this_copy = *this;

	auto setDeltas = [&](const Operator& left, const Operator& right, bool sc_type) {
		Momentum copy_momentum = left.momentum;
		if (sc_type) copy_momentum.flipMomentum();

		for (size_t j = 1; j < left.indizes.size(); j++)
		{
			this->delta_indizes.push_back(std::make_pair(left.indizes[j], right.indizes[j]));
			this_copy.delta_indizes.push_back(this->delta_indizes.back());
		}

		this->delta_momenta.push_back(std::make_pair(copy_momentum, right.momentum));
		copy_momentum.add_Q = !(copy_momentum.add_Q);
		this_copy.delta_momenta.push_back(std::make_pair(copy_momentum, right.momentum));
	};

	for (size_t i = 0; i < temporary_operators.size(); i += 2)
	{
		if (LEFT.isDaggered == RIGHT.isDaggered) {
			if (L_SPIN == R_SPIN) return false;
			if (LEFT.isDaggered) { // c^+ c^+
				if (L_SPIN == DOWN) {
					std::cerr << "c^+ c^+: Left spin is down while right isn't. Did you forget to sort the terms?" << std::endl;
					throw;
				}
				this->delta_indizes.push_back(std::make_pair(L_SPIN, UP));
				this_copy.delta_indizes.push_back(this->delta_indizes.back());
				this->delta_indizes.push_back(std::make_pair(R_SPIN, DOWN));
				this_copy.delta_indizes.push_back(this->delta_indizes.back());
				// Due to the dagger we need to swap left and right
				setDeltas(RIGHT, LEFT, true);

				this->operators.push_back(WickOperator("f", LEFT.momentum));
				if (LEFT.indizes.size() > 1) {
					this->operators.back().indizes = std::vector<std::string>(LEFT.indizes.begin() + 1, LEFT.indizes.end());
				}
				this_copy.operators.push_back(this->operators.back());
				this_copy.operators.back().type = "\\eta";
			}
			else { // cc
				if (L_SPIN == UP) {
					std::cerr << "c c: Left spin is up while right isn't. Did you forget to sort the terms?" << std::endl;
					throw;
				}
				this->delta_indizes.push_back(std::make_pair(L_SPIN, DOWN));
				this_copy.delta_indizes.push_back(this->delta_indizes.back());
				this->delta_indizes.push_back(std::make_pair(R_SPIN, UP));
				this_copy.delta_indizes.push_back(this->delta_indizes.back());
				setDeltas(LEFT, RIGHT, true);

				this->operators.push_back(WickOperator("f", RIGHT.momentum));
				if (RIGHT.indizes.size() > 1) {
					this->operators.back().indizes = std::vector<std::string>(RIGHT.indizes.begin() + 1, RIGHT.indizes.end());
				}
				this_copy.operators.push_back(this->operators.back());
				this_copy.operators.back().type = "\\eta";
			}
		}
		else {
			// c^+ c
			if (L_SPIN == UP && R_SPIN == DOWN) return false;
			if (L_SPIN == DOWN && R_SPIN == UP) return false;

			this->delta_indizes.push_back(std::make_pair(L_SPIN, R_SPIN));
			this_copy.delta_indizes.push_back(this->delta_indizes.back());
			// Left and right are swapped due to the definition of g
			setDeltas(RIGHT, LEFT, false);

			this->operators.push_back(WickOperator("n", LEFT.momentum, LEFT.indizes));
			this_copy.operators.push_back(this->operators.back());
			this_copy.operators.back().type = "g";
		}
	}
	reciever.push_back(this_copy);
	return true;
}

bool WickTerm::setDeltas()
{
	for (auto& delta : delta_momenta)
	{
		for (auto it = delta.first.momentum_list.begin(); it != delta.first.momentum_list.end(); )
		{
			int index = delta.second.isUsed(it->second);
			if (index < 0) {
				++it;
				continue;
			}

			int remainder = delta.second.momentum_list[index].first - it->first;
			if (remainder == 0) {
				delta.second.momentum_list.erase(delta.second.momentum_list.begin() + index);
				it = delta.first.momentum_list.erase(it);
				continue;
			}

			delta.second.momentum_list[index].first = remainder;
			it = delta.first.momentum_list.erase(it);
		}
		if (delta.first.add_Q == delta.second.add_Q) {
			delta.first.add_Q = false;
			delta.second.add_Q = false;
		}
		if (delta.first.momentum_list.size() == 0) {
			if (delta.second.momentum_list.size() == 0) continue;
			std::swap(delta.first, delta.second);
		}
		if (delta.first.add_Q) {
			delta.first.add_Q = false;
			delta.second.add_Q = !(delta.second.add_Q);
		}
		if (delta.first.momentum_list.front().first < 0) {
			delta.first.flipMomentum();
			delta.second.flipMomentum();
		}
		if (delta.first.momentum_list.size() > 1 && delta.second.momentum_list.size() == 0) {
			delta.second.momentum_list.push_back(delta.first.momentum_list[1]);
			delta.first.momentum_list.erase(delta.first.momentum_list.begin() + 1);
		}
	}

	for (auto it = delta_momenta.begin(); it != delta_momenta.end(); )
	{
		if (it->first.momentum_list.size() == 0 && it->second.momentum_list.size() == 0) {
			// 0 = Q can never be achieved
			if (it->first.add_Q != it->second.add_Q) return false;
			it = delta_momenta.erase(it);
		}
		else {
			++it;
		}
	}

	// Set all deltas up to the same notation
	for (auto& delta : delta_momenta) {
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
			}
			if (!foundCandidate) index = 0;

			if (delta.second.momentum_list[index].first > 0) {
				delta.second.flipMomentum();
			}
			delta.first.momentum_list.push_back(delta.second.momentum_list[index]);
			delta.first.flipMomentum();
			if (abs(delta.first.momentum_list[0].first) != 1) std::cerr << "Not yet implemented! " << delta.first << std::endl;
			delta.second.momentum_list.erase(delta.second.momentum_list.begin() + index);
		}

		if (abs(delta.first.momentum_list[0].first) != 1) std::cerr << "Not yet implemented! " << delta.first << std::endl;
		for (auto& op : operators) {
			op.momentum.replaceOccurances(delta.first.momentum_list[0].second, delta.second);
		}
		for (auto& coeff : coefficients) {
			coeff.momentum.replaceOccurances(delta.first.momentum_list[0].second, delta.second);
		}
	}
	for (auto& delta : delta_indizes) {
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
	for (size_t i = 0; i < delta_momenta.size(); i++)
	{
	setDeltas_outerLoop_momentum:

		for (size_t j = i + 1; j < delta_momenta.size(); j++)
		{
			if (pair_equal_allow_permutation(delta_momenta[i], delta_momenta[j])) {
				delta_momenta.erase(delta_momenta.begin() + j);
				goto setDeltas_outerLoop_momentum;
			}

			auto delta_buffer = delta_momenta[j];
			delta_buffer.first.flipMomentum();
			delta_buffer.second.flipMomentum();
			if (pair_equal_allow_permutation(delta_momenta[i], delta_buffer)) {
				delta_momenta.erase(delta_momenta.begin() + j);
				goto setDeltas_outerLoop_momentum;
			}
		}
	}
	for (size_t i = 0; i < delta_indizes.size(); i++)
	{
	setDeltas_outerLoop_index:

		for (size_t j = i + 1; j < delta_indizes.size(); j++)
		{
			if (pair_equal_allow_permutation(delta_indizes[i], delta_indizes[j])) {
				delta_indizes.erase(delta_indizes.begin() + j);
				goto setDeltas_outerLoop_index;
			}
		}
	}

	// Erase delta_k,k etc
	for (auto it = delta_momenta.begin(); it != delta_momenta.end();)
	{
		// k = k + Q can never be achieved
		if (it->first.differsOnlyInQ(it->second)) return false;

		if (it->first == it->second) {
			it = delta_momenta.erase(it);
		}
		else {
			++it;
		}
	}
	for (auto it = delta_indizes.begin(); it != delta_indizes.end();)
	{
		// UP can never be DOWN and vice versa
		if (it->first == UP && it->second == DOWN) return false;
		if (it->first == DOWN && it->second == UP) return false;

		if (it->first == it->second) {
			it = delta_indizes.erase(it);
		}
		else {
			++it;
		}
	}

	return true;
}

void WickTerm::computeSums()
{
	auto changeAllIndizes = [&](const std::string& replaceWhat, const std::string& replaceWith) {
		for (auto& op : operators) {
			for (auto it = op.indizes.begin(); it != op.indizes.end(); ++it)
			{
				if (*it == replaceWhat) {
					*it = replaceWith;
				}
			}
		}
		for (auto& coeff : coefficients) {
			for (auto it = coeff.indizes.begin(); it != coeff.indizes.end(); ++it)
			{
				if (*it == replaceWhat) {
					*it = replaceWith;
				}
			}
		}
	};

	for (size_t i = 0; i < sum_indizes.size(); i++)
	{
	computeSums_index_loop:

		for (size_t j = 0; j < delta_indizes.size(); j++)
		{
			if (delta_indizes[j].first == sum_indizes[i]) {
				changeAllIndizes(sum_indizes[i], delta_indizes[j].second);
				sum_indizes.erase(sum_indizes.begin() + i);
				delta_indizes.erase(delta_indizes.begin() + j);
				goto computeSums_index_loop;
			}
			else if (delta_indizes[j].second == sum_indizes[i]) {
				changeAllIndizes(sum_indizes[i], delta_indizes[j].first);
				sum_indizes.erase(sum_indizes.begin() + i);
				delta_indizes.erase(delta_indizes.begin() + j);
				goto computeSums_index_loop;
			}
		}
	}

	auto changeAllMomenta = [&](const char replaceWhat, const Momentum& replaceWith) {
		for (auto& op : operators) {
			op.momentum.replaceOccurances(replaceWhat, replaceWith);
		}
		for (auto& coeff : coefficients) {
			coeff.momentum.replaceOccurances(replaceWhat, replaceWith);
		}
	};

	for (size_t i = 0; i < sum_momenta.size(); i++)
	{
	computeSums_momentum_loop:

		for (size_t j = 0; j < delta_momenta.size(); j++)
		{
			if (delta_momenta[j].first.momentum_list[0].second == sum_momenta[i]) {
				changeAllMomenta(sum_momenta[i], delta_momenta[j].second);
				if (abs(delta_momenta[j].first.momentum_list[0].first) != 1) std::cerr << "Not yet implemented! " << delta_momenta[j].first << std::endl;

				sum_momenta.erase(sum_momenta.begin() + i);
				delta_momenta.erase(delta_momenta.begin() + j);
				goto computeSums_momentum_loop;
			}
			else {
				int index = delta_momenta[j].second.isUsed(sum_momenta[i]);
				if (index < 0) continue;

				Momentum buffer(delta_momenta[j].second.momentum_list[index].second, delta_momenta[j].second.momentum_list[index].first);
				if (abs(buffer.momentum_list[0].first) != 1) std::cerr << "Not yet implemented! " << buffer << std::endl;
				delta_momenta[j].second.momentum_list.erase(delta_momenta[j].second.momentum_list.begin() + index);
				delta_momenta[j].second -= delta_momenta[j].first;

				if (buffer.momentum_list[0].first > 0) {
					delta_momenta[j].second.flipMomentum();
				}
				changeAllMomenta(sum_momenta[i], delta_momenta[j].second);

				sum_momenta.erase(sum_momenta.begin() + i);
				delta_momenta.erase(delta_momenta.begin() + j);
				goto computeSums_momentum_loop;
			}
		}
	}
}

void WickTerm::discardZeroMomenta()
{
	for (auto& op : operators) {
		for (auto it = op.momentum.momentum_list.begin(); it != op.momentum.momentum_list.end();) {
			if (it->first == 0) {
				it = op.momentum.momentum_list.erase(it);
			}
			else {
				++it;
			}
		}
	}
	for (auto& coeff : coefficients) {
		for (auto it = coeff.momentum.momentum_list.begin(); it != coeff.momentum.momentum_list.end();) {
			if (it->first == 0) {
				it = coeff.momentum.momentum_list.erase(it);
			}
			else {
				++it;
			}
		}
	}
}

void WickTerm::renameSums()
{
	const char name_list[3] = { 'q', 'p', 'r' };
	const char buffer_list[3] = { 'x', 'y', 'z' };
	for (size_t i = 0; i < sum_momenta.size(); i++)
	{
		if (i >= 3) {
			std::cerr << "More than 3 momenta, time to implement this..." << std::endl;
			break;
		}
		if (sum_momenta[i] == name_list[i]) continue;

		for (auto& op : operators) {
			op.momentum.replaceOccurances(sum_momenta[i], Momentum(buffer_list[i]));
		}
		for (auto& coeff : coefficients) {
			coeff.momentum.replaceOccurances(sum_momenta[i], Momentum(buffer_list[i]));
		}
		sum_momenta[i] = name_list[i];
	}

	for (size_t i = 0; i < sum_momenta.size(); i++)
	{
		for (auto& op : operators) {
			op.momentum.replaceOccurances(buffer_list[i], Momentum(name_list[i]));
		}
		for (auto& coeff : coefficients) {
			coeff.momentum.replaceOccurances(buffer_list[i], Momentum(name_list[i]));
		}
	}

	for (const auto& sum : sum_momenta)
	{
		for (auto& op : operators) {
			int index = op.momentum.isUsed(sum);
			if (index < 0) continue;
			if (op.momentum.momentum_list.size() == 1) break;

			Momentum buffer = op.momentum;
			buffer.flipMomentum();
			buffer.momentum_list[index].first *= -1;
			buffer.momentum_list[index].second = 'x';

			for (auto& op2 : operators) {
				op2.momentum.replaceOccurances(sum, buffer);
				op2.momentum.replaceOccurances('x', Momentum(sum));
			}
			for (auto& coeff : coefficients) {
				coeff.momentum.replaceOccurances(sum, buffer);
				coeff.momentum.replaceOccurances('x', Momentum(sum));
			}
		}
	}
	this->discardZeroMomenta();
}

void WickTerm::sort()
{
	std::map<std::string, int> sort_map;
	sort_map["n"] = 0;
	sort_map["g"] = 1;
	sort_map["f"] = 2;
	sort_map["\\eta"] = 3;

	for (size_t i = 0; i < operators.size(); i++)
	{
		for (size_t j = i + 1; j < operators.size(); j++)
		{
			if (sort_map.find(operators[i].type)->second > sort_map.find(operators[j].type)->second) {
				std::swap(operators[i], operators[j]);
			}
		}
	}
}

void cleanWicks(std::vector<WickTerm>& terms)
{
	for (auto it = terms.begin(); it != terms.end();) {
		if (!(it->setDeltas())) {
			it = terms.erase(it);
			continue;
		}
		it->discardZeroMomenta();
		it->computeSums();

		if (!(it->setDeltas())) {
			it = terms.erase(it);
			continue;
		}
		it->discardZeroMomenta();
		it->renameSums();
		it->sort();
		++it;
	}

	// remove duplicates
	TODO!

	// removes any terms that have a 0 prefactor
	for (auto it = terms.begin(); it != terms.end();)
	{
		if (it->multiplicity == 0) {
			it = terms.erase(it);
		}
		else {
			++it;
		}
	}
}

std::ostream& operator<<(std::ostream& os, const WickOperator& op)
{
	os << "\\langle " << op.type << "_{ " << op.momentum << ", ";
	for (const auto& index : op.indizes) {
		os << index << " ";
	}
	os << "} \\rangle";
	return os;
}
std::ostream& operator<<(std::ostream& os, const std::vector<WickOperator>& ops)
{
	for (const auto& op : ops) {
		os << op << " ";
	}
	return os;
}
std::ostream& operator<<(std::ostream& os, const WickTerm& term)
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
	os << term.coefficients << " ";
	for (const auto& delta : term.delta_momenta) {
		os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
	}
	for (const auto& delta : term.delta_indizes) {
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
std::ostream& operator<<(std::ostream& os, const std::vector<WickTerm>& terms)
{
	for (std::vector<WickTerm>::const_iterator it = terms.begin(); it != terms.end(); ++it)
	{
		os << "\t&" << *it;
		if (it != terms.end() - 1) {
			os << " \\\\";
		}
		os << "\n";
	}
	return os;
}

WickOperator::WickOperator(const std::string& _type, const Momentum& _momentum, const std::vector<std::string>& _indizes)
	: type(_type), momentum(_momentum), indizes(_indizes) {}
WickOperator::WickOperator(const std::string& _type, const Momentum& _momentum, const std::string& _index)
	: type(_type), momentum(_momentum), indizes(1, _index) {}
WickOperator::WickOperator()
	: type(""), momentum(), indizes() {}