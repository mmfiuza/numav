// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "hash-functions.hpp"

namespace std
{
    uint64_t hash<pair<uint64_t, uint64_t>>::operator()(
        const pair<uint64_t, uint64_t> key
    ) const {
        return (get<0UL>(key) << (4UL * sizeof(uint64_t))) + get<1UL>(key);
    }

    uint64_t hash<tuple<uint64_t, uint64_t>>::operator()(
        const tuple<uint64_t, uint64_t> key
    ) const {
        return (get<0UL>(key) << (4UL * sizeof(uint64_t))) + get<1UL>(key);
    }

    uint64_t hash<std::vector<uint64_t>>::operator()(
        const std::vector<uint64_t> vec
    ) const {
        uint64_t seed = vec.size();
        for (auto& i : vec) {
            seed ^= i + 0x9e3779b9 + (seed << 6UL) + (seed >> 2UL);
        }
        return seed;
    }
}
