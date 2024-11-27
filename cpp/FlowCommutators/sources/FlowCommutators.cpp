#include <SymbolicOperators/Term.hpp>

using namespace SymbolicOperators;

int main(int argc, char** argv) {
    const Term H_C(IntFractional(1, 2), Coefficient("V", Momentum('q')),
			SumContainer{ MomentumSum({ 'k', 'l', 'q' }), IndexSum({ Index::Sigma, Index::SigmaPrime }) },
			std::vector<Operator>({
				Operator('k', 1, false, Index::Sigma, true),
				Operator('l', 1, false, Index::SigmaPrime, true),
				Operator(momentum_pairs({ std::make_pair(1, 'k'), std::make_pair(-1, 'q') }), Index::SigmaPrime, false),
				Operator(momentum_pairs({ std::make_pair(1, 'l'), std::make_pair(1, 'q') }), Index::Sigma, false),
				}));
    
    const Term CUT_eta_annihilation(-1, Coefficient("A", MomentumList({ 'P', 'Q' })), 
                SumContainer{ MomentumSum({ 'P', 'Q' }), IndexSum( Index::Sigma ) },
                std::vector<Operator>({
                    Operator::Boson(Momentum('Q'), false),
                    Operator(momentum_pairs({ std::make_pair(1, 'P'), std::make_pair(1, 'Q') }), Index::Sigma, true),
                    Operator('P', 1, false, Index::Sigma, false)
                }));
    const Term CUT_eta_creation(1, Coefficient("B", MomentumList({ 'P', 'Q' })), 
                SumContainer{ MomentumSum({ 'P', 'Q' }), IndexSum( Index::Sigma ) },
                std::vector<Operator>({
                    Operator::Boson(Momentum('Q', -1), false),
                    Operator(momentum_pairs({ std::make_pair(1, 'P'), std::make_pair(1, 'Q') }), Index::Sigma, true),
                    Operator('P', 1, false, Index::Sigma, false)
                }));
    const std::vector<Term> CUT_eta { CUT_eta_annihilation, CUT_eta_creation };
    std::cout << "\\begin{align*}\n\t \\eta = " 
              << CUT_eta 
              << "\\end{align*}" << std::endl;
    
    std::vector<Term> commutation_result;
    commutator(commutation_result, CUT_eta, H_C);
    cleanUp(commutation_result);
    for(auto& term : commutation_result){
        term.renameMomenta('r', 'k');
        term.renameMomenta('s', 'l');
    }
    std::cout << "\\begin{align*}\n\t [\\eta, H_\\mathrm{C}] = " 
              << commutation_result 
              << "\\end{align*}" << std::endl;
    return 0;
}