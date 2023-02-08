#include "Term.hpp"

void Term::print() const{
    std::cout << multiplicity << " \\cdot ";
    if(!sum_indizes.empty()){
        std::cout << "\\sum_{ ";
        for(const auto& index : sum_indizes){
            std::cout << index << " ";
        }
        std::cout << "}";
    }
    if(!sum_momenta.empty()){
        std::cout << "\\sum_{ ";
        for(const auto& momentum : sum_momenta){
            std::cout << momentum << " ";
        }
        std::cout << "}";
    }
    std::cout << this->coefficient << " \\cdot ";
    for(const auto& delta : delta_momentum){
        std::cout << "\\delta_{" << delta.first << ", " << delta.second << "} ";
    }
    for(const auto& delta : delta_index){
        std::cout << "\\delta_{" << delta.first << ", " << delta.second << "} ";
    }

    if(isIdentity()){
        std::cout << " \\mathbb{1} " << std::endl;
        return;
    }
    for(const auto& op : operators){
        std::cout << op << " ";
    }
    std::cout << std::endl;
}

void normalOrder(std::vector<Term>& terms){
    for(size_t t = 0; t < terms.size();){
        size_t n = terms[t].operators.size();
        size_t new_n;
        while(n > 1){
            n = 0;
            for (size_t i = 1; i < terms[t].operators.size(); i++)
            {
                if(!(terms[t].operators[i - 1].isDaggered) && (terms[t].operators[i].isDaggered)){
                    new_n = i;
                    // Swap cc^+
                    terms[t].multiplicity *= -1;
                    std::swap(terms[t].operators[i - 1], terms[t].operators[i]);

                    // Add a new term where cc^+ is replaced by the appropriate delta
                    Term new_term(terms[t]);
                    if(new_term.operators[i - 1].indizes.size() != new_term.operators[i].indizes.size()) {
                        std::cerr << "Operators do not have the same index count." << std::endl;
                        throw;
                    }
                    for (size_t c = 0; c < new_term.operators[i - 1].indizes.size(); c++)
                    {
                        // if the indizes are not the same we emplace a delta
                        // otherwise no action is required
                        if(new_term.operators[i - 1].indizes[c] != new_term.operators[i - 1].indizes[c]){
                            new_term.delta_index.push_back(
                                std::make_pair(new_term.operators[i - 1].indizes[c], new_term.operators[i].indizes[c]));
                        }
                    }
                    if(new_term.operators[i-1].momentum != new_term.operators[i].momentum){
                        new_term.delta_momentum.push_back(
                            std::make_pair(new_term.operators[i-1].momentum, new_term.operators[i].momentum)
                        );
                    }

                    new_term.operators.erase(new_term.operators.begin() + i - 1, new_term.operators.begin() + i);
                } 
                else if(terms[t].operators[i - 1] == terms[t].operators[i]){
                    // two identical fermion operators = 0
                    terms.erase(terms.begin() + t);
                    goto outerLoop;
                }
            }
            n = new_n;
        }

        ++t;
        outerLoop:
        
    } 
}