#pragma once
#include <vector>
#include <string>

#include "Operator.hpp"

struct Coefficient{
    std::string name;
    Momentum momentum;
    // Contains all indizes, standard: first index = spin, all others arbitrary, e.g. orbitals, bands etc
    std::vector<std::string> indizes;
    bool isDaggered;

    Coefficient() : name(""), momentum(), indizes(), isDaggered(false) {};
    Coefficient(std::string _name, const Momentum& _momentum, const std::vector<std::string>& _indizes, bool _isDaggered)
        : name(_name), momentum(_momentum), indizes(_indizes), isDaggered(_isDaggered){};
    Coefficient(std::string _name, char _momentum, bool add_Q, const std::vector<std::string>& _indizes, bool _isDaggered)
        : name(_name), momentum(_momentum, add_Q), indizes(_indizes), isDaggered(_isDaggered){ };
};

std::ostream& operator<<(std::ostream& os, const Coefficient& coeff){
    os << coeff.name;
    if(!coeff.indizes.empty()){
        os << "_{ ";
        for(const auto& index : coeff.indizes){
            os << index << " ";
        }
        os << "}";
    }
    if(coeff.isDaggered){
        os << "^*";
    }
    if(!coeff.momentum.momentum_list.empty()){
        os << " ( " << coeff.momentum << " )";
    }

    return os;
};

class Term{
private:
    int multiplicity;
    Coefficient coefficient;
    std::vector<char> sum_momenta;
    std::vector<std::string> sum_indizes;
    std::vector<Operator> operators;
    // symbolises the Kronecker delta
    std::vector<std::pair<Momentum, Momentum>> delta_momentum;
    std::vector<std::pair<std::string, std::string>> delta_index;

public:
    Term(int _multiplicity, Coefficient _coefficient, const std::vector<char>& _sum_momenta, const std::vector<std::string>& _sum_indizes, const std::vector<Operator>& _operators=std::vector<Operator>())
        : multiplicity(_multiplicity), coefficient(_coefficient), sum_momenta(_sum_momenta), sum_indizes(_sum_indizes), operators(_operators) {};

    Term(int _multiplicity, Coefficient _coefficient, const std::vector<char>& _sum_momenta, const std::vector<Operator>& _operators=std::vector<Operator>())
        : multiplicity(_multiplicity), coefficient(_coefficient), sum_momenta(_sum_momenta), operators(_operators) {};

    Term(int _multiplicity, Coefficient _coefficient, const std::vector<std::string>& _sum_indizes, const std::vector<Operator>& _operators=std::vector<Operator>())
        : multiplicity(_multiplicity), coefficient(_coefficient), sum_indizes(_sum_indizes), operators(_operators) {};

    Term(int _multiplicity, Coefficient _coefficient, const std::vector<Operator>& _operators=std::vector<Operator>())
        : multiplicity(_multiplicity), coefficient(_coefficient), operators(_operators) {};

    inline bool isIdentity() const {
        return this->operators.empty();
    }

    void print() const;
    friend void normalOrder(std::vector<Term>& terms);
};