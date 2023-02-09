#define append_vector(a, b) a.insert(a.end(), b.begin(), b.end())
#include "Operator.hpp"

Momentum& Momentum::operator+=(const Momentum& rhs)
{
	this->add_Q = (rhs.add_Q != this->add_Q);
	append_vector(this->momentum_list, rhs.momentum_list);

	std::vector<std::pair<int, char>>::iterator it = momentum_list.begin();
	std::vector<std::pair<int, char>>::iterator jt = momentum_list.begin();
	while (it != momentum_list.end()) {
		jt = it + 1;
		while (jt != momentum_list.end()) {
			if (it->second == jt->second) {
				it->first += jt->first;
				jt = this->momentum_list.erase(jt);
			}
			else {
				++jt;
			}
		}
		if (it != momentum_list.end()) {
			++it;
		}
	}
	return *this;
}

Momentum& Momentum::operator-=(const Momentum& rhs)
{
	this->add_Q = (rhs.add_Q != this->add_Q);
	const size_t old_size = this->momentum_list.size();
	append_vector(this->momentum_list, rhs.momentum_list);

	for (size_t i = old_size; i < this->momentum_list.size(); i++)
	{
		this->momentum_list[i].first *= -1;
	}
	std::vector<std::pair<int, char>>::iterator it = momentum_list.begin();
	std::vector<std::pair<int, char>>::iterator jt = momentum_list.begin();
	while (it != momentum_list.end()) {
		jt = it + 1;
		while (jt != momentum_list.end()) {
			if (it->second == jt->second) {
				it->first += jt->first;
				jt = this->momentum_list.erase(jt);
			}
			else {
				++jt;
			}
		}
		if (it != momentum_list.end()) {
			++it;
		}
	}
	return *this;
}

void Momentum::addInPlace(const Momentum& rhs)
{
	(*this) += rhs;
}

void Momentum::replaceOccurances(const char replaceWhat, const Momentum& replaceWith)
{
	for (size_t i = 0; i < momentum_list.size(); ++i) {
		if (momentum_list[i].second == replaceWhat) {
			auto buffer = replaceWith;
			buffer.multiplyMomentum(momentum_list[i].first);
			this->momentum_list.erase(momentum_list.begin() + i);

			this->addInPlace(buffer);
		}
	}
}

std::ostream& operator<<(std::ostream& os, const Momentum& momentum)
{
	if (momentum.momentum_list.empty()) {
		os << "0";
		return os;
	}
	for (momentum_pairs::const_iterator it = momentum.momentum_list.begin(); it != momentum.momentum_list.end(); ++it)
	{
		if (it != momentum.momentum_list.begin() && it->first > 0) {
			os << " +";
		}
		else if (it->first == -1) {
			os << " -";
		}
		else if (abs(it->first) != 1) {
			os << it->first;
		}
		os << it->second;
	}

	return os;
}

std::ostream& operator<<(std::ostream& os, const Operator& op)
{
	os << "c_{ " << op.momentum << ", ";
	for (const auto& index : op.indizes) {
		os << index << " ";
	}
	os << "}";
	if (op.isDaggered) {
		os << "^\\dagger ";
	}
	return os;
}