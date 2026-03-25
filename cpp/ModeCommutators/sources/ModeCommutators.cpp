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

std::vector<Term> get_hamiltonian()
{
    const Term H_Kin(1, Coefficient("\\epsilon_0", Momentum('q')), SumContainer{ MomentumSum({ 'q' }), Index::Sigma },
		std::vector<Operator>({
			Operator('q', 1, false, Index::Sigma, true), Operator('q', 1, false, Index::Sigma, false)
			})
        );
    const Term H_Ph(1, Coefficient::RealInversionSymmetric("g", MomentumList({ 'q', 'p' }), std::function<void(Coefficient&)>([](Coefficient& coeff){ coeff.momenta.sort(); })),
		SumContainer{ MomentumSum({ 'p', 'q' }) },
		std::vector<Operator>({
			c_k_dagger.with_momentum('q'), c_minus_k_dagger.with_momentum('q'),
			c_minus_k.with_momentum('p'), c_k.with_momentum('p') })
		);
    return { H_Kin, H_Ph};
}

int main(int argc, char** argv) {
	std::vector<Term> hamiltonian = get_hamiltonian();
	
	const std::vector<Term> pc_term = std::vector<Term>({
					Term(1, std::vector<Operator>({ c_k_dagger, c_minus_k_dagger }))
			});
	const std::vector<Term> higgs_term = std::vector<Term>({
					Term(1, std::vector<Operator>({ c_minus_k, c_k })),
					Term(1, std::vector<Operator>({ c_k_dagger, c_minus_k_dagger }))
			});
	const std::vector<Term> number_term = std::vector<Term>({
					Term(1, std::vector<Operator>({ c_minus_k_dagger, c_minus_k })),
					Term(1, std::vector<Operator>({ c_k_dagger, c_k }))
			});

	const std::vector<Term> delta_k_operator = std::vector<Term>({
		Term(1, Coefficient::RealInversionSymmetric("g", MomentumList({ 'k', 't' }), std::function<void(Coefficient&)>([](Coefficient& coeff){ coeff.momenta.sort(); })),
			SumContainer{ MomentumSum({ 't' }) },
			std::vector<Operator>({
				c_k_dagger.with_momentum('t'), c_minus_k_dagger.with_momentum('t')
		}))
	});


	std::vector<Term> higgs_commutator = commutator(hamiltonian, higgs_term);
	clean_up(higgs_commutator);
	std::ranges::sort(higgs_commutator, [](const Term& l, const Term& r) {
		if (l.operators.size() < r.operators.size()) return true;
		if (l.operators.size() > r.operators.size()) return false;
		if (!l.operators.front().is_daggered && r.operators.front().is_daggered) return true;
		return false;
		});
	std::cout << "\\begin{align*}\n\t[H, c_{-k \\downarrow } c_{k \\uparrow } + c_{k \\uparrow}^\\dagger c_{-k \\downarrow}^\\dagger ] = "
		<< higgs_commutator
		<< "\\end{align*}" << std::endl;

	std::vector<Term> number_commutator = commutator(hamiltonian, number_term);
	clean_up(number_commutator);
	std::ranges::sort(number_commutator, [](const Term& l, const Term& r) {
		if (l.operators.size() < r.operators.size()) return true;
		if (l.operators.size() > r.operators.size()) return false;
		if (!l.operators.front().is_daggered && r.operators.front().is_daggered) return true;
		return false;
		});
	std::cout << "\\begin{align*}\n\t[H, c_{k \\uparrow}^\\dagger c_{k \\uparrow } + c_{-k \\downarrow}^\\dagger c_{-k \\downarrow} ] = "
		<< number_commutator
		<< "\\end{align*}" << std::endl;
	
	
	std::vector<Term> delta_commutator = commutator(hamiltonian, delta_k_operator);
	clean_up(delta_commutator);
	std::ranges::sort(delta_commutator, [](const Term& l, const Term& r) {
		if (l.operators.size() < r.operators.size()) return true;
		if (l.operators.size() > r.operators.size()) return false;
		if (!l.operators.front().is_daggered && r.operators.front().is_daggered) return true;
		return false;
		});
	std::cout << "\\begin{align*}\n\t[H, \\hat{\\Delta}_{k}^\\dagger ] = "
		<< delta_commutator
		<< "\\end{align*}" << std::endl;
	
	std::vector<Term> pc_commutator = commutator(hamiltonian, pc_term);
	clean_up(pc_commutator);
	std::ranges::sort(pc_commutator, [](const Term& l, const Term& r) {
		if (l.operators.size() < r.operators.size()) return true;
		if (l.operators.size() > r.operators.size()) return false;
		if (!l.operators.front().is_daggered && r.operators.front().is_daggered) return true;
		return false;
		});
	std::cout << "\\begin{align*}\n\t[H, c_{k \\uparrow}^\\dagger c_{-k \\downarrow}^\\dagger ] = "
		<< pc_commutator
		<< "\\end{align*}" << std::endl;

	return 0;
}