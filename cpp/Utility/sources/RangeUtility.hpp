#pragma once
#include <iterator>

namespace Utility{
    template<class Vector>
    void duplicate_n_inplace(Vector& target, size_t n){
        const size_t original_size = target.size();
        target.reserve(original_size * n);
        for(size_t i = 0U; i < n; ++i){
            std::copy_n(target.begin(), original_size, std::back_inserter(target));
        }
    }

    template<class Vector>
    void append_vector(Vector& target, const Vector& source){
        target.insert(target.end(), source.begin(), source.end());
    }
}