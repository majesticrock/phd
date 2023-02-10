#define append_vector(a, b) a.insert(a.end(), b.begin(), b.end())
#include "Operator.hpp"

Momentum& Momentum::operator+=(const Momentum& rhs)
{
	this->add_Q = (rhs.add_Q != this->add_Q);
	bool foundOne = false;
	for (size_t i = 0; i < rhs.momentum_list.size(); i++)
	{
		foundOne = false;
		for (size_t j = 0; j < this->momentum_list.size(); j++)
		{
			if(rhs.momentum_list[i].second == this->momentum_list[j].second){
				foundOne = true;
				this->momentum_list[j].first += rhs.momentum_list[i].first;
			}
		}
		if(!foundOne){
			this->momentum_list.push_back(rhs.momentum_list[i]);
		}
	}
	
	return *this;
}

Momentum& Momentum::operator-=(const Momentum& rhs)
{
	this->add_Q = (rhs.add_Q != this->add_Q);
	bool foundOne = false;
	for (size_t i = 0; i < rhs.momentum_list.size(); i++)
	{
		foundOne = false;
		for (size_t j = 0; j < this->momentum_list.size(); j++)
		{
			if(rhs.momentum_list[i].second == this->momentum_list[j].second){
				foundOne = true;
				this->momentum_list[j].first -= rhs.momentum_list[i].first;
			}
		}
		if(!foundOne){
			this->momentum_list.push_back(std::make_pair(-rhs.momentum_list[i].first, rhs.momentum_list[i].second));
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

			(*this) +=buffer;
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