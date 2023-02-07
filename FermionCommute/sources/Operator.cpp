#include "Operator.hpp"

void Momentum::addInPlace(const Momentum& rhs)
{
    this->add_Q = (rhs.add_Q != this->add_Q);
    this->momentum_list.insert(std::end(this->momentum_list), std::begin(rhs.momentum_list), std::end(rhs.momentum_list));
    
    std::vector<std::pair<int, char>>::iterator it = momentum_list.begin();
    std::vector<std::pair<int, char>>::iterator jt = momentum_list.begin();
    while(it != momentum_list.end()){
        jt = it + 1;
        while(jt != momentum_list.end()){
            if(it->second == jt->second){
                it->first += jt->first;
                jt = this->momentum_list.erase(jt);
            }else{
                ++jt;
            }
        }
        if(it != momentum_list.end()){
            ++it;
        }
    }
}