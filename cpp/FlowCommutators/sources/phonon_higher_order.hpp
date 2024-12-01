#pragma once
#include <SymbolicOperators/Term.hpp>
#include <vector>
#include <utility>

inline void phonon_higher_order() {
	using namespace SymbolicOperators;
	const std::vector<Term> eta_prime({
		Term(1, Coefficient("A_2", MomentumList({ 'K', 'L', 'P', 'Q' }), IndexWrapper{}, false, false),
			SumContainer{ MomentumSum({ 'K', 'L', 'P', 'Q' }), IndexSum({ Index::Sigma, Index::SigmaPrime }) },
			std::vector<Operator>({
				Operator::Boson(Momentum('P', -1), true),
				Operator(momentum_pairs({ std::make_pair(1, 'K'), std::make_pair(1, 'P')}), Index::Sigma, true),
				Operator('L', 1, false, Index::SigmaPrime, true),
				Operator(momentum_pairs({ std::make_pair(1, 'L'), std::make_pair(-1, 'Q') }), Index::SigmaPrime, false),
				Operator(momentum_pairs({ std::make_pair(1, 'K'), std::make_pair(1, 'Q') }), Index::Sigma, false),
			})),
		Term(1, Coefficient("B_2", MomentumList({ 'K', 'L', 'P', 'Q' }), IndexWrapper{}, false, false),
			SumContainer{ MomentumSum({ 'K', 'L', 'P', 'Q' }), IndexSum({ Index::Sigma, Index::SigmaPrime }) },
			std::vector<Operator>({
				Operator::Boson(Momentum('P', 1), false),
				Operator(momentum_pairs({ std::make_pair(1, 'K'), std::make_pair(1, 'P')}), Index::Sigma, true),
				Operator('L', 1, false, Index::SigmaPrime, true),
				Operator(momentum_pairs({ std::make_pair(1, 'L'), std::make_pair(-1, 'Q') }), Index::SigmaPrime, false),
				Operator(momentum_pairs({ std::make_pair(1, 'K'), std::make_pair(1, 'Q') }), Index::Sigma, false),
			}))
		});
	const std::vector<Term> H_original({
		// c^+ c
		Term(1, Coefficient("\\epsilon", Momentum('k')),
			SumContainer{ MomentumSum({'k'}), IndexSum({ Index::Sigma }) },
			std::vector<Operator>({
				Operator(Momentum('k'), Index::Sigma, true),
				Operator(Momentum('k'), Index::Sigma, false),
			})),
		// b^+ b
		Term(1, Coefficient("\\omega", Momentum('q')),
			SumContainer{ MomentumSum({'q'}), {} },
			std::vector<Operator>({
				Operator::Boson(Momentum('q'), true),
				Operator::Boson(Momentum('q'), false)
			})),
		// el-ph interaction
		Term(1, Coefficient("M", MomentumList({ 'k', 'q' }), IndexWrapper{}, false, false),
			SumContainer{ MomentumSum({ 'k', 'q' }), IndexSum(Index::Sigma) },
			std::vector<Operator>({
				Operator::Boson(Momentum('q', -1), true),
				Operator(momentum_pairs({ std::make_pair(1, 'k'), std::make_pair(1, 'q') }), Index::Sigma, true),
				Operator('k', 1, false, Index::Sigma, false)
			})),
		Term(-1,Coefficient("M", MomentumList({ 'k', 'q' }), IndexWrapper{}, false, false, true),
			SumContainer{ MomentumSum({ 'k', 'q' }), IndexSum(Index::Sigma) },
			std::vector<Operator>({
				Operator::Boson(Momentum('q'), false),
				Operator(momentum_pairs({ std::make_pair(1, 'k'), std::make_pair(1, 'q') }), Index::Sigma, true),
				Operator('k', 1, false, Index::Sigma, false)
			})),
		// el-el interaction
		Term(1, Coefficient("V_ph", MomentumList({ 'k', 'l', 'q' }), IndexWrapper{}, false, false),
			SumContainer{ MomentumSum({ 'k', 'l', 'q' }), IndexSum({ Index::Sigma, Index::SigmaPrime }) },
			std::vector<Operator>({
				Operator('k', 1, false, Index::Sigma, true),
				Operator('l', 1, false, Index::SigmaPrime, true),
				Operator(momentum_pairs({ std::make_pair(1, 'l'), std::make_pair(-1, 'q') }), Index::SigmaPrime, false),
				Operator(momentum_pairs({ std::make_pair(1, 'k'), std::make_pair(1, 'q') }), Index::Sigma, false),
			}))
		});

	std::vector<Term> commutation_result;
	commutator(commutation_result, eta_prime, H_original);
	cleanUp(commutation_result);

	std::erase_if(commutation_result, [](const Term& term) -> bool {
		if (term.operators.size() > 5U) return true;
		return false;
		});
	std::sort(commutation_result.begin(), commutation_result.end(), [](const Term& l, const Term& r) {
		return l.operators.size() < r.operators.size();
		});
	

	std::cout << "\\begin{align*}\n\t[\\eta', H] = "
		<< commutation_result
		<< ". \\end{align*}" << std::endl;
}