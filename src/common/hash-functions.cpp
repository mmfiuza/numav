// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "hash-functions.hpp"

namespace std
{
    // define a hash function for std::pair<size_t,size_t>
    size_t hash<pair<size_t, size_t>>::operator()(
        const pair<size_t, size_t>& key
    ) const {
        return (get<0>(key) << (4 * sizeof(size_t))) + get<1>(key);
    }

    // define a hash function for std::tuple<size_t,size_t>
    size_t hash<tuple<size_t, size_t>>::operator()(
        const tuple<size_t, size_t>& key
    ) const {
        return (get<0>(key) << (4 * sizeof(size_t))) + get<1>(key);
    }
}
