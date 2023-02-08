#include "Term.hpp"

int main(int argc, char** argv){
    Momentum k('k');
    Operator c_k(k, std::vector<std::string>(), false);
    Operator c_k_dagger(k, std::vector<std::string>(), true);

    Term term(1, Coefficient(), std::vector<Operator>({ 
        c_k_dagger, c_k
    }));

    term.print();
    return 0;
}