// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

// define a hash function for std::pair<size_t,size_t>
namespace std {
    template<>
    struct hash<std::pair<size_t,size_t>> {
        size_t operator()(const std::pair<size_t,size_t>& key) const {
            return (std::get<0>(key) << 4*sizeof(size_t)) + std::get<1>(key);
        }
    };
}

// define a hash function for std::tuple<size_t,size_t>
namespace std {
    template<>
    struct hash<std::tuple<size_t,size_t>> {
        size_t operator()(const std::tuple<size_t,size_t>& key) const {
            return (std::get<0>(key) << 4*sizeof(size_t)) + std::get<1>(key);
        }
    };
}
