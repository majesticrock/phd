#pragma once
#include <vector>
#include <string>
#include <iostream>

struct Momentum{
        // total momentum is then:    sum_i pair_i.first * pair_i.second
        std::vector<std::pair<int, char>> momentum_list;
        bool add_Q;

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
        : indizes(_indizes), isDaggered(_isDaggered){
            this->momentum = {std::make_pair<int, char>(1, _momentum), add_Q };
        };
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
                foundOne == true;
            }
        }
        if(!foundOne) return false;
    }
    return true;
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