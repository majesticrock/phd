#include <SymbolicOperators/Term.hpp>
#include <vector>

using namespace SymbolicOperators;

int main(int argc, char** argv) {
	const Term H_C(IntFractional(1, 4), 
		Coefficient("V", MomentumList({ 'k', 'l', 'q' }), IndexWrapper{}, false, false),
		SumContainer{ MomentumSum({ 'k', 'l', 'q' }), IndexSum({ Index::Sigma, Index::SigmaPrime }) },
		std::vector<Operator>({
			Operator('k', 1, false, Index::Sigma, true),
			Operator('l', 1, false, Index::SigmaPrime, true),
			Operator(momentum_pairs({ std::make_pair(1, 'l'), std::make_pair(-1, 'q') }), Index::SigmaPrime, false),
			Operator(momentum_pairs({ std::make_pair(1, 'k'), std::make_pair(1, 'q') }), Index::Sigma, false),
		}));
	const Term H_C_mixed(IntFractional(1, 4), 
		Coefficient("V", MomentumList({ Momentum('l'), Momentum('k'), Momentum('q', -1) }), IndexWrapper{}, false, false),
		SumContainer{ MomentumSum({ 'k', 'l', 'q' }), IndexSum({ Index::Sigma, Index::SigmaPrime }) },
		std::vector<Operator>({
			Operator('k', 1, false, Index::Sigma, true),
			Operator('l', 1, false, Index::SigmaPrime, true),
			Operator(momentum_pairs({ std::make_pair(1, 'l'), std::make_pair(-1, 'q') }), Index::SigmaPrime, false),
			Operator(momentum_pairs({ std::make_pair(1, 'k'), std::make_pair(1, 'q') }), Index::Sigma, false),
		}));
	const std::vector<Term> H_sym{H_C, H_C_mixed};
	
	const Term CUT_eta_creation(1, 
		Coefficient("A", MomentumList({ 'P', 'Q' }), IndexWrapper{}, false, false),
		SumContainer{ MomentumSum({ 'P', 'Q' }), IndexSum(Index::GeneralSpin_S) },
		std::vector<Operator>({
			Operator::Boson(Momentum('Q', -1), true),
			Operator(momentum_pairs({ std::make_pair(1, 'P'), std::make_pair(1, 'Q') }), Index::GeneralSpin_S, true),
			Operator('P', 1, false, Index::GeneralSpin_S, false)
		}));
	const Term CUT_eta_annihilation(-1, 
		Coefficient("B", MomentumList({ 'P', 'Q' }), IndexWrapper{}, false, false), 
		SumContainer{ MomentumSum({ 'P', 'Q' }), IndexSum( Index::GeneralSpin_S ) },
		std::vector<Operator>({
			 Operator::Boson(Momentum('Q'), false),
			 Operator(momentum_pairs({ std::make_pair(1, 'P'), std::make_pair(1, 'Q') }), Index::GeneralSpin_S, true),
			 Operator('P', 1, false, Index::GeneralSpin_S, false)
		 }));
	
	const std::vector<Term> CUT_eta { CUT_eta_creation, CUT_eta_annihilation };
	std::cout << "The Coulomb interaction term is" << std::endl;
	std::cout << "\\begin{align*}\n\tH_\\mathrm{C} = "
		<< H_sym
		<< ",\n\\end{align*}" << std::endl;
	std::cout << "Here, we used the symmetry of the potential $V(k, l, q) = V(l, k, q)$ and $V(k, l, q) = V(k, l, -q)$."
		<< "We can even define\n\\begin{equation}\n\tV^\\mathrm{sym} (k, l, q) \\equiv (1/2) \\cdot (V(k, l, q) + V(l, k, -q))\n\\end{equation}\n"
		<< "because this quantity will appear throughout the flow." << std::endl;
	std::cout << "Our generator for the CUT is " << std::endl;
	std::cout << "\\begin{align*}\n\t\\eta = " 
			  << CUT_eta 
			  << ", \\end{align*}" << std::endl;
	std::cout << "where\n\\begin{align}\n\tA(P,Q) &= M_{P,Q} \\sgn ( \\epsilon_{P+Q} - \\epsilon_P + \\omega_{-Q} ) \\\\\n"
		<< "\tB(P,Q) &= M_{P+Q,-Q}^* \\sgn ( \\epsilon_{P+Q} - \\epsilon_P - \\omega_Q )\n\\end{align}\n" << std::endl;

	std::vector<Term> commutation_result;
	commutator(commutation_result, CUT_eta, H_sym);
	cleanUp(commutation_result);
	for (auto& term : commutation_result) {
		term.renameMomenta('q', 'x');
		term.renameMomenta('p', 'k');
		term.renameMomenta('r', 'l');
		term.renameMomenta('s', 'q');
		term.renameMomenta('x', 'p');

		// Third terms of the
		if (term.operators[3].momentum.momentum_list.size() == 3U) {
			Momentum& base = term.operators[3].momentum;
			Momentum replacement = Momentum('x') + Momentum('p');
			term.transform_momentum_sum(base.momentum_list[0].second, replacement, 'x');
			term.renameMomenta('x', 'l');
		}
		else if (term.operators[4].momentum.momentum_list.size() == 3U) {
			Momentum& base = term.operators[4].momentum;
			Momentum replacement = Momentum('x') + Momentum('p');
			term.transform_momentum_sum(base.momentum_list[0].second, replacement, 'x');
			term.renameMomenta('x', 'k');
		}
		// We may have terms like c_k^+ c_l+p^+ ... instead of c_k+p^+ c_l^+ ...
		// Those are actually identical by swapping k <-> l, q <-> -q, and sigma <-> sigma'
		// And lastly swapping the two creation operators, and swapping the two annihilation operators.
		// We perfom this action now
		if (term.operators[2].momentum.momentum_list.size() == 2) {
			// Starting at 1, because the boson does not have a spin
			for (size_t i = 1U; i < term.operators.size(); ++i) {
				if (term.operators[i].first_index() == Index::Sigma) {
					term.operators[i].set_first_index(Index::SigmaPrime);
				}
				else if (term.operators[i].first_index() == Index::SigmaPrime) {
					term.operators[i].set_first_index(Index::Sigma);
				}
				else {
					throw std::runtime_error("Expected only Sigma and SigmaPrime as indizes!");
				}
			}
			term.renameMomenta('k', 'x');
			term.renameMomenta('l', 'k');
			term.renameMomenta('x', 'l');
			term.invert_momentum_sum('q');
			std::swap(term.operators[1], term.operators[2]);
			std::swap(term.operators[3], term.operators[4]);

			// Using the symmetry of the Coulomb potential V(q) = V(-q)
			// The second coefficient is the Coulomb one
			//Coefficient& coulomb = term.coefficients[1];
			//assert(coulomb.name == "V");
			//assert(coulomb.momenta.back().momentum_list.size() == 1);
			//coulomb.momenta.back().momentum_list[0].first = std::abs(coulomb.momenta.back().momentum_list[0].first);
		}
		// Sort the sum indizes to be identical
		std::sort(term.sums.momenta.begin(), term.sums.momenta.end());
	}
	clear_duplicates(commutation_result);

	std::cout << "The first commutation yields" << std::endl;
	std::cout << "\\begin{align*}\n\t[\\eta, H_\\mathrm{C}] = " 
			  << commutation_result 
			  << ". \\end{align*}" << std::endl;
	std::cout << "Here, we recognize our $V^\\mathrm{sym}$! Symplifying the expression to" << std::endl;
	for (auto& term : commutation_result) {
		Coefficient& coeff = term.coefficients[1];
		if (coeff.momenta[0].uses('l')) {
			std::swap(coeff.momenta[0], coeff.momenta[1]);
		}
		if (coeff.momenta[2] == Momentum('q', -1)) {
			coeff.momenta[2].flipMomentum();
		}
		coeff.name = "V^\\mathrm{sym}";
	}
	clear_duplicates(commutation_result);
	std::cout << "\\begin{align*}\n\t[\\eta, H_\\mathrm{C}] = " 
			  << commutation_result 
			  << ", \\end{align*}" << std::endl;
	std::cout << "which has no two-electron interaction terms." 
		<< " Thus, we include the resulting terms into our Hamiltonian as " << std::endl;

	const std::vector<Term> H_prime(
		{
			Term(1, Coefficient("C_1", MomentumList({ 'k', 'l', 'p', 'q' }), IndexWrapper{}, false, false),
			SumContainer{ MomentumSum({ 'k', 'l', 'p', 'q' }), IndexSum({ Index::Sigma, Index::SigmaPrime }) },
			std::vector<Operator>({
				Operator::Boson(Momentum('p', -1), true),
				Operator(momentum_pairs({ std::make_pair(1, 'k'), std::make_pair(1, 'p')}), Index::Sigma, true),
				Operator('l', 1, false, Index::SigmaPrime, true),
				Operator(momentum_pairs({ std::make_pair(1, 'l'), std::make_pair(-1, 'q') }), Index::SigmaPrime, false),
				Operator(momentum_pairs({ std::make_pair(1, 'k'), std::make_pair(1, 'q') }), Index::Sigma, false),
			})),
			Term(1, Coefficient("C_2", MomentumList({ 'k', 'l', 'p', 'q' }), IndexWrapper{}, false, false),
			SumContainer{ MomentumSum({ 'k', 'l', 'p', 'q' }), IndexSum({ Index::Sigma, Index::SigmaPrime }) },
			std::vector<Operator>({
				Operator::Boson(Momentum('p', 1), false),
				Operator(momentum_pairs({ std::make_pair(1, 'k'), std::make_pair(1, 'p')}), Index::Sigma, true),
				Operator('l', 1, false, Index::SigmaPrime, true),
				Operator(momentum_pairs({ std::make_pair(1, 'l'), std::make_pair(-1, 'q') }), Index::SigmaPrime, false),
				Operator(momentum_pairs({ std::make_pair(1, 'k'), std::make_pair(1, 'q') }), Index::Sigma, false),
			})),
		});


	std::cout << "\\begin{align*}\n\tH'= "
		<< H_prime
		<< ", \\end{align*}" << std::endl;

	std::vector<Term> second_commutation_result;
	commutator(second_commutation_result, CUT_eta, H_prime);
	cleanUp(second_commutation_result);
	for (auto& term : second_commutation_result) {
		term.renameMomenta('p', 'l');
		term.renameMomenta('q', 'p');
		term.renameMomenta('r', 'k');
		term.renameMomenta('s', 'q');
		term.renameMomenta('t', 'Q');
	}

	std::cout << "We commute again with $H'$ and obtain" << std::endl;
	std::cout << "\\begin{align*}\n\t[\\eta, H'] = "
		<< second_commutation_result
		<< "\\end{align*}" << std::endl;
	std::cout << "Here, we want to omit all terms with bosonic creation/annihilation operators, claiming them to be small. " 
		<< "We proceed in the same manner with hexatic terms." << std::endl;
	std::cout << "It is also nice to directly sort the summations so that the operators appear as they do in the Coulomb interaction. This yields" << std::endl;

	std::erase_if(second_commutation_result, [](const Term& term) { return (term.contains_boson() || term.operators.size() == 6U); });

	for (auto& term : second_commutation_result) {
		{
			Momentum& first_momentum = term.operators[0].momentum;
			if (first_momentum != Momentum('k')) {
				if (first_momentum.momentum_list.size() == 1U) {
					const char original = first_momentum.momentum_list[0].second;
					term.renameMomenta(original, 'x');
					term.renameMomenta('k', original);
					term.renameMomenta('x', 'k');
				} 
				else {
					int k_pos = first_momentum.isUsed('k');
					assert(k_pos > -1);
					Momentum replacement = Momentum('x') - first_momentum + Momentum(first_momentum.momentum_list[k_pos]);
					term.transform_momentum_sum(first_momentum.momentum_list[k_pos].second, replacement, 'x');
					term.renameMomenta('x', 'k');
				}
			}
		}
		{
			Momentum& second_momentum = term.operators[1].momentum;
			if (second_momentum != Momentum('l')) {
				int l_pos = second_momentum.isUsed('l');
				assert(l_pos > -1);
				Momentum replacement = Momentum('x') - second_momentum + Momentum(second_momentum.momentum_list[l_pos]);
				term.transform_momentum_sum(second_momentum.momentum_list[l_pos].second, replacement, 'x');
				term.renameMomenta('x', 'l');
			}
		}
		{
			Momentum& third_momentum = term.operators[2].momentum;
			int q_pos = third_momentum.isUsed('q');
			assert(q_pos > -1);
			int l_pos = third_momentum.isUsed('l');
			assert(l_pos > -1);
			Momentum replacement = Momentum('x') - third_momentum + Momentum(third_momentum.momentum_list[q_pos]) + Momentum(third_momentum.momentum_list[l_pos]);
			if(third_momentum.momentum_list[q_pos].first  < 0) {
				replacement.flipMomentum();
			}
			term.transform_momentum_sum(third_momentum.momentum_list[q_pos].second, replacement, 'x');
			term.renameMomenta('x', 'q');

			q_pos = third_momentum.isUsed('q');
			assert(q_pos == 1);
			if (third_momentum.momentum_list[q_pos].first > 0) {
				term.invert_momentum_sum('q');
			}
		}
	}
	{
		// We are allowed to change the p in any manner we like - to obtain the same notation across all terms we do the following:
		Term* term_ptr = &(second_commutation_result[0]);
		// \sum_{ p l k q } \sum_{ \sigma \sigma' } A( k+q, -k+p )  C_2( p, l, k-p, q )   c_{ k, \sigma }^\dagger  c_{ l, \sigma' }^\dagger  c_{ l-q, \sigma' } c_{ k+q, \sigma }
		Momentum replacement = Momentum('x') + Momentum('k');
		term_ptr->transform_momentum_sum('p', replacement, 'x');
		term_ptr->renameMomenta('x', 'p');

		term_ptr = &(second_commutation_result[1]);
		// \sum_{ p l k q } \sum_{ \sigma \sigma' } A( l-q, -k+p )  C_2( p, l, k-p, k-p+q )   c_{ k, \sigma }^\dagger  c_{ l, \sigma' }^\dagger  c_{ l-q, \sigma' } c_{ k+q, \sigma }
		term_ptr->transform_momentum_sum('p', replacement, 'x');
		term_ptr->renameMomenta('x', 'p');

		term_ptr = &(second_commutation_result[3]);
		// \sum_{ p l k q } \sum_{ \sigma \sigma' } B( k+l-p, -k+p )  C_1( p, k+l-p, k-p, k-p+q )   c_{ k, \sigma }^\dagger  c_{ l, \sigma' }^\dagger  c_{ l-q, \sigma' } c_{ k+q, \sigma } 
		term_ptr->transform_momentum_sum('p', replacement, 'x');
		term_ptr->renameMomenta('x', 'p');
	}

	std::cout << "\\begin{align*}\n\t[\\eta, H'] \\approx "
		<< second_commutation_result
		<< "\\end{align*}" << std::endl;

	std::cout << "Inserting the definitions for $A$ and $B$ yields" << std::endl;
	std::cout << "\\begin{align*}\n\t\\partial_{\\lambda} V(\\lambda; k, l, q) \\approx - \\sum_{p} \\bigg\\{ ";
	
	std::array<Coefficient, 4> M_lambdas;
	M_lambdas.fill(Coefficient("M_{\\lambda}"));
	std::array<std::array<Coefficient, 3>, 4> in_signums;
	in_signums.fill({ Coefficient("\\epsilon"), Coefficient("\\epsilon"), Coefficient("\\omega") });
	for (size_t i = 0U; i < second_commutation_result.size(); ++i) {
		// The first coefficient is always A or B
		const Coefficient& term_coeff = second_commutation_result[i].coefficients[0];
		assert(term_coeff.name == "A" || term_coeff.name == "B");
		M_lambdas[i].translational_invariance = false;
		M_lambdas[i].momenta.resize(2);
		if (term_coeff.name == "A") {
			M_lambdas[i].momenta[0] = term_coeff.momenta[0];
			M_lambdas[i].momenta[1] = term_coeff.momenta[1];
		} 
		else {
			M_lambdas[i].momenta[0] = term_coeff.momenta[0] + term_coeff.momenta[1];
			M_lambdas[i].momenta[1] = -term_coeff.momenta[1];
			M_lambdas[i].is_daggered = true;
		}

		in_signums[i][0].momenta.push_back(term_coeff.momenta[0] + term_coeff.momenta[1]);
		in_signums[i][1].momenta.push_back(term_coeff.momenta[0]);
		in_signums[i][2].momenta.push_back(term_coeff.name == "A" ? -term_coeff.momenta[1] : term_coeff.momenta[1]);
		
		if(i > 0) {
			std::cout << "\t+ ";
		}
		std::cout << "&" << M_lambdas[i] << " \\sgn \\left[ ";
		for(int j = 0; j < 3; ++j) {
			if (j == 1) {
				std::cout << " - ";
			}
			else if (j == 2) {
				std::cout << (term_coeff.name == "A" ? "+" : "-");
			}
			std::cout << in_signums[i][j];
		}
		std::cout << " \\right]" << second_commutation_result[i].coefficients[1];
		if(i < second_commutation_result.size() - 1) std::cout << "\\\\\n";
	}
	std::cout << " \\bigg\\}\n\\end{align*}" << std::endl;
	return 0;
}