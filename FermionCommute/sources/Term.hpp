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

    Coefficient(std::string name, const Momentum& _momentum, const std::vector<std::string>& _indizes, bool _isDaggered)
        : name(_name), momentum(_momentum), indizes(_indizes), isDaggered(_isDaggered){};
    Coefficient(std::string name, char _momentum, bool add_Q, const std::vector<std::string>& _indizes, bool _isDaggered)
        : name(_name), indizes(_indizes), isDaggered(_isDaggered){
        this->momentum = {std::make_pair<int, char>(1, _momentum), add_Q };
     };
};

class Term{
private:
    int multiplicity;
    Coefficient coefficient;
    std::vector<char> sum_momenta;
    std::vector<std::string> sum_indizes;
    std::vector<Operator> operators;

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

    friend void normalOrder(std::vector<Term>& terms);
};