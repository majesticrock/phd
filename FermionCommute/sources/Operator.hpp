#pragma once
#include <vector>
#include <string>
#include <iostream>

typedef std::vector<std::pair<int, char>> momentum_pairs ;
struct Momentum{
    // total momentum is then:    sum_i pair_i.first * pair_i.second
    momentum_pairs momentum_list;
    bool add_Q;

    Momentum() : momentum_list(), add_Q(false) {};
    explicit Momentum(const char value, bool Q=false) : momentum_list(1, std::make_pair(1, value)), add_Q(Q) {};
    Momentum(momentum_pairs& _momenta, bool Q) : momentum_list(_momenta), add_Q(Q){};
    void addInPlace(const Momentum& rhs);
};

struct Operator{
    Momentum momentum;
    // Contains all indizes, standard: first index = spin, all others arbitrary, e.g. orbitals, bands etc
    std::vector<std::string> indizes;
    bool isDaggered;

    Operator(const Momentum& _momentum, const std::vector<std::string>& _indizes, bool _isDaggered)
        : momentum(_momentum), indizes(_indizes), isDaggered(_isDaggered){};
    Operator(char _momentum, bool add_Q, const std::vector<std::string>& _indizes, bool _isDaggered)
        : momentum(_momentum, add_Q), indizes(_indizes), isDaggered(_isDaggered){};
};

std::ostream& operator<<(std::ostream& os, const Momentum& momentum){
    if(momentum.momentum_list.empty()){
        os << "0";
        return os;
    }
    for (momentum_pairs::const_iterator it = momentum.momentum_list.begin(); it != momentum.momentum_list.end(); ++it)
    {
        if(it != momentum.momentum_list.begin() && it->first > 0){
            os << " +";
        } else if(it->first == -1){
            os << " -";
        }
        else if(abs(it->first) != 1){
            os << it->first;
        }
        os << it->second;
    }
    
    return os;
};

std::ostream& operator<<(std::ostream& os, const Operator& op){
    os << "c_{ " << op.momentum << ", ";
    for(const auto& index : op.indizes){
        os << index << " ";
    }
    os << "}";
    if(op.isDaggered){ 
        os << "^\\dagger ";
    }
    return os;
};

inline bool operator==(const Momentum& lhs, const Momentum& rhs){
    if(lhs.add_Q != rhs.add_Q) return false;
    if(lhs.momentum_list.size() != rhs.momentum_list.size()) return false;
    bool foundOne = false;
    for (size_t i = 0; i < lhs.momentum_list.size(); i++)
    {
        foundOne = false;
        for (size_t j = 0; j < rhs.momentum_list.size(); j++)
        {
            if(lhs.momentum_list[i] == rhs.momentum_list[j]){
                foundOne = true;
            }
        }
        if(!foundOne) return false;
    }
    return true;
}

inline bool operator!=(const Momentum& lhs, const Momentum& rhs){
    return !(lhs == rhs);
}

inline bool operator==(const Operator& lhs, const Operator& rhs){
    if(lhs.isDaggered != rhs.isDaggered){
        return false;
    }
    if(lhs.indizes.size() != rhs.indizes.size()){
        // Should never occur, but better to check anyways
        std::cout << "Warning: Number of indizes between two operators does not match!" << std::endl;
        return false;
    }
    for (size_t i = 0; i < lhs.indizes.size(); i++)
    {
        if(lhs.indizes[i] != rhs.indizes[i]){
            return false;
        }
    }
    return (lhs.momentum == rhs.momentum);
}

inline bool operator!=(const Operator& lhs, const Operator& rhs){
    return !(lhs == rhs);
}