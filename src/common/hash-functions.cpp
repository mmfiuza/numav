// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "hash-functions.hpp"

namespace std
{
    size_t hash<pair<size_t, size_t>>::operator()(
        const pair<size_t, size_t>& key
    ) const {
        return (get<0UL>(key) << (4UL * sizeof(size_t))) + get<1UL>(key);
    }

    size_t hash<tuple<size_t, size_t>>::operator()(
        const tuple<size_t, size_t>& key
    ) const {
        return (get<0UL>(key) << (4UL * sizeof(size_t))) + get<1UL>(key);
    }

    size_t hash<std::vector<size_t>>::operator()(
        const std::vector<size_t>& vec
    ) const {
        size_t seed = vec.size();
        for (auto& i : vec) {
            seed ^= i + 0x9e3779b9 + (seed << 6UL) + (seed >> 2UL);
        }
        return seed;
    }
}
