#pragma once
#include <vector>
#include <iostream>

typedef std::vector<std::pair<int, char>> momentum_pairs;
struct Momentum {
	// total momentum is then:    sum_i pair_i.first * pair_i.second
	momentum_pairs momentum_list;
	bool add_Q;

	Momentum();
	explicit Momentum(const char value, int plus_minus = 1, bool Q = false);
	explicit Momentum(const momentum_pairs& _momenta, bool Q = false);

	Momentum& operator+=(const Momentum& rhs);
	Momentum& operator-=(const Momentum& rhs);

	inline void multiplyMomentum(int factor) {
		for (auto& m : momentum_list) {
			m.first *= factor;
		}
	};
	inline void flipMomentum() {
		multiplyMomentum(-1);
	};
	inline int isUsed(const char value) const {
		for (int i = 0; i < momentum_list.size(); ++i) {
			if (momentum_list[i].second == value) return i;
		}
		return -1;
	};

	void addInPlace(const Momentum& rhs);
	void replaceOccurances(const char replaceWhat, const Momentum& replaceWith);
};

inline Momentum operator+(Momentum lhs, const Momentum& rhs) {
	lhs += rhs;
	return lhs;
}
inline Momentum operator-(Momentum lhs, const Momentum& rhs) {
	lhs -= rhs;
	return lhs;
}
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

std::ostream& operator<<(std::ostream& os, const Momentum& momentum);