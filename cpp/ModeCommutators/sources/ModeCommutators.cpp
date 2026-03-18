#include <mrock/symbolic_operators/Momentum.hpp>
#include <mrock/symbolic_operators/Operator.hpp>
#include <mrock/symbolic_operators/Term.hpp>

#include <vector>

using namespace mrock::symbolic_operators;

const Momentum base_k = Momentum('k');
const Operator c_k = Operator{ base_k, Index::SpinUp, false };
const Operator c_minus_k = Operator{ -base_k, Index::SpinDown, false };
const Operator c_k_dagger = Operator{ base_k, Index::SpinUp, true };
const Operator c_minus_k_dagger = Operator{ -base_k, Index::SpinDown, true };

std::vector<Term> hamiltonian(char momentum_1, char momentum_2)
{
    const Term H_Kin(1, Coefficient("\\epsilon_0", Momentum(momentum_1)), SumContainer{ MomentumSum({ momentum_1 }), Index::Sigma },
		std::vector<Operator>({
			Operator(momentum_1, 1, false, Index::Sigma, true), Operator(momentum_1, 1, false, Index::Sigma, false)
			})
        );
    const Term H_Ph(-1, Coefficient("g"),
		SumContainer{ MomentumSum({ momentum_2, momentum_1 }) },
		std::vector<Operator>({
			c_k_dagger.with_momentum(momentum_1), c_minus_k_dagger.with_momentum(momentum_1),
			c_minus_k.with_momentum(momentum_2), c_k.with_momentum(momentum_2) })
		);
    return { H_Kin, H_Ph};
}

int main(int argc, char** argv) {
	std::vector<Term> hamiltonian1 = hamiltonian('q', 'p');
	std::vector<Term> hamiltonian2 = hamiltonian('r', 's');
	
	std::vector<Term> phase_term = std::vector<Term>({
					Term(1, std::vector<Operator>({ c_minus_k, c_k })),
					Term(-1, std::vector<Operator>({ c_k_dagger, c_minus_k_dagger }))
			});

	std::vector<Term> inner_commutator = commutator(hamiltonian1, phase_term);
	clean_up(inner_commutator);
	
	std::vector<Term> commutation_result = commutator(hamiltonian2, inner_commutator);
	clean_up(commutation_result);

	std::cout << "\\begin{align*}\n\tH = "
		<< hamiltonian1
		<< "\\end{align*}" << std::endl;
	
	std::ranges::sort(commutation_result, [](const Term& l, const Term& r) {
		if (l.operators.size() < r.operators.size()) return true;
		if (l.operators.size() > r.operators.size()) return false;
		if (!l.operators.front().is_daggered && r.operators.front().is_daggered) return true;
		return false;
		});

	std::cout << "\\begin{align*}\n\t[H, [H, c_{-k \\downarrow } c_{k \\uparrow } - c_{k \\uparrow}^\\dagger c_{-k \\downarrow}^\\dagger ]] = "
		<< commutation_result
		<< "\\end{align*}" << std::endl;
	return 0;
}