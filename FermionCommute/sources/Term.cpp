#include "Term.hpp"

void normalOrder(std::vector<Term>& terms){
    for(std::vector<Term>::iterator term = terms.begin(); term != terms.end();){
        size_t n = term->operators.size();
        size_t new_n;
        while(n > 1){
            n = 0;
            for (size_t i = 1; i < term->operators.size(); i++)
            {
                if(!(term->operators[i - 1].isDaggered) && (term->operators[i].isDaggered)){
                    new_n = i;
                    // Swap cc^+
                    term->multiplicity *= -1;
                    std::swap(term->operators[i - 1], term->operators[i]);

                    // Add a new term where cc^+ is replaced by the appropriate delta
                }
            }
            n = new_n;
        }

        ++term;
    } 
}