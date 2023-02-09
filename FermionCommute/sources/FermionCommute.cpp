#include "Term.hpp"

int main(int argc, char** argv) {
	Momentum k('k');
	Momentum minus_k('k', -1);
	Momentum l('l');
	Momentum minus_l('l', -1);
	Momentum diff(momentum_pairs({ std::make_pair(1, 'p'), std::make_pair(-1, 'k') }), false);

	Operator c_k(k, std::vector<std::string>({ UP }), false);
	Operator c_k_dagger(k, std::vector<std::string>({ UP }), true);
	Operator c_l(l, std::vector<std::string>({ UP }), false);
	Operator c_l_dagger(l, std::vector<std::string>({ UP }), true);

	Operator c_minus_k(minus_k, std::vector<std::string>({ DOWN }), false);
	Operator c_minus_k_dagger(minus_k, std::vector<std::string>({ DOWN }), true);
	Operator c_minus_l(minus_l, std::vector<std::string>({ DOWN }), false);
	Operator c_minus_l_dagger(minus_l, std::vector<std::string>({ DOWN }), true);

	Operator c_diff(diff, std::vector<std::string>({ UP }), false);
	Operator c_diff_dagger(diff, std::vector<std::string>({ UP }), true);

	Term left(1, Coefficient(), std::vector<Operator>({
		c_minus_k, c_k
		}));
	Term right(1, Coefficient(), std::vector<Operator>({
		c_l_dagger, c_minus_l_dagger
		}));

	std::vector<Term> terms;
	commutator(terms, left, right);
	cleanUp(terms);

	std::cout << "\\begin{align*}\n" << "[" << left.toStringWithoutPrefactor() << ", " << right.toStringWithoutPrefactor() << "] =" 
		<< terms << "\\end{align*}" << std::endl;

	return 0;
}