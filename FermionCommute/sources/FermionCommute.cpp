#include "Term.hpp"

int main(int argc, char** argv) {
	Operator c_k('k', 1, false, UP, false);
	Operator c_minus_k('k', -1, false, DOWN, false);

	Operator c_l_dagger('l', 1, false, UP, true);
	Operator c_minus_l_dagger('l', -1, false, DOWN, true);

	Term H_T(1, Coefficient("\\epsilon_0", 'q'), std::vector<char>({ 'q' }), std::vector<std::string>({ "\\sigma" }), std::vector<Operator>({
		Operator('q', 1, false, "\\sigma", true), Operator('q', 1, false, "\\sigma", false)
		}));

	Term H_U(1, Coefficient("\\frac{U}{N}"), std::vector<char>({ 'r', 'p', 'q' }), std::vector<Operator>({
		Operator('r', 1, false, UP, true), Operator('p', 1, false, DOWN, true),
		Operator(momentum_pairs({ std::make_pair(1, 'p'), std::make_pair(-1, 'q') }), DOWN, false),
		Operator(momentum_pairs({ std::make_pair(1, 'r'), std::make_pair(1, 'q') }), UP, false),
		}));

	std::vector<Term> H = { H_T, H_U };

	Term right(1, Coefficient(), std::vector<Operator>({
		c_minus_k, c_k
		}));

	Term left(1, Coefficient(), std::vector<Operator>({
		c_l_dagger, c_minus_l_dagger
		}));

	std::vector<Term> commute_with_H, terms;
	commutator(commute_with_H, H, right);
	cleanUp(commute_with_H);

	std::cout << "\\begin{align*}\n\t[ H, " << right.toStringWithoutPrefactor() << "] = " << commute_with_H << "\\end{align*}" << std::endl;

	std::vector<WickTerm> wicks;
	for (const auto& term : commute_with_H) {
		term.wick(wicks);
	}
	cleanWicks(wicks);

	std::cout << "\\begin{align*}\n\t\\langle[ H, " << right.toStringWithoutPrefactor() << "] \\rangle = "
		<< wicks << "\\end{align*}" << std::endl;
	return 0;

	commutator(terms, left, commute_with_H);
	cleanUp(terms);

	std::cout << "\\begin{align*}\n\t[" << left.toStringWithoutPrefactor() << ", [ H, " << right.toStringWithoutPrefactor() << "] ] = "
		<< terms << "\\end{align*}" << std::endl;

	return 0;
}