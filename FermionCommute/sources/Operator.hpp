#pragma once
#include <vector>
#include <string>
#include <iostream>

const std::string UP = "\\uparrow";
const std::string DOWN = "\\downarrow";

typedef std::vector<std::pair<int, char>> momentum_pairs;
struct Momentum {
	// total momentum is then:    sum_i pair_i.first * pair_i.second
	momentum_pairs momentum_list;
	bool add_Q;

	Momentum() : momentum_list(), add_Q(false) {};
	explicit Momentum(const char value, int plus_minus = 1, bool Q = false) : momentum_list(1, std::make_pair(plus_minus, value)), add_Q(Q) {};
	explicit Momentum(const momentum_pairs& _momenta, bool Q = false) : momentum_list(_momenta), add_Q(Q) {};

	Momentum& operator+=(const Momentum& rhs);
	Momentum& operator-=(const Momentum& rhs);

	inline void multiplyMomentum(int factor) {
		for (auto& m : momentum_list) {
			m.first *= factor;
		}
	}
	inline void flipMomentum() {
		multiplyMomentum(-1);
	};
	inline int isUsed(const char value) const {
		for (int i = 0; i < momentum_list.size(); ++i) {
			if (momentum_list[i].second == value) return i;
		}
		return -1;
	}

	void addInPlace(const Momentum& rhs);
	void replaceOccurances(const char replaceWhat, const Momentum& replaceWith);
};

struct Operator {
	Momentum momentum;
	// Contains all indizes, standard: first index = spin, all others arbitrary, e.g. orbitals, bands etc
	std::vector<std::string> indizes;
	bool isDaggered;

	Operator(const Momentum& _momentum, const std::vector<std::string>& _indizes, bool _isDaggered)
		: momentum(_momentum), indizes(_indizes), isDaggered(_isDaggered) {};
	Operator(const Momentum& _momentum, const std::string& index, bool _isDaggered)
		: momentum(_momentum), indizes(1, index), isDaggered(_isDaggered) {};
	Operator(const momentum_pairs& _momentum, const std::string& index, bool _isDaggered)
		: momentum(_momentum), indizes(1, index), isDaggered(_isDaggered) {};
	Operator(char _momentum, bool add_Q, const std::vector<std::string>& _indizes, bool _isDaggered)
		: momentum(_momentum, add_Q), indizes(_indizes), isDaggered(_isDaggered) {};
	Operator(char _momentum, int sign, bool add_Q, const std::string& index, bool _isDaggered)
		: momentum(_momentum, sign, add_Q), indizes(1, index), isDaggered(_isDaggered) {};
};

std::ostream& operator<<(std::ostream& os, const Momentum& momentum);

std::ostream& operator<<(std::ostream& os, const Operator& op);

inline Momentum operator+(Momentum lhs, const Momentum& rhs) {
	lhs += rhs;
	return lhs;
};
inline Momentum operator-(Momentum lhs, const Momentum& rhs) {
	lhs -= rhs;
	return lhs;
};

inline bool operator==(const Momentum& lhs, const Momentum& rhs) {
	if (lhs.add_Q != rhs.add_Q) return false;
	if (lhs.momentum_list.size() != rhs.momentum_list.size()) return false;
	bool foundOne = true;
	for (size_t i = 0; i < lhs.momentum_list.size(); i++)
	{
		foundOne = false;
		for (size_t j = 0; j < rhs.momentum_list.size(); j++)
		{
			if (lhs.momentum_list[i] == rhs.momentum_list[j])
				foundOne = true;
		}
		if (!foundOne) return false;
	}
	return true;
}

inline bool operator!=(const Momentum& lhs, const Momentum& rhs) {
	return !(lhs == rhs);
}

inline bool operator==(const Operator& lhs, const Operator& rhs) {
	if (lhs.isDaggered != rhs.isDaggered) {
		return false;
	}
	if (lhs.indizes.size() != rhs.indizes.size()) {
		// Should never occur, but better to check anyways
		std::cout << "Warning: Number of indizes between two operators does not match!" << std::endl;
		return false;
	}
	for (size_t i = 0; i < lhs.indizes.size(); i++)
	{
		if (lhs.indizes[i] != rhs.indizes[i]) {
			return false;
		}
	}
	return (lhs.momentum == rhs.momentum);
}

inline bool operator!=(const Operator& lhs, const Operator& rhs) {
	return !(lhs == rhs);
}