// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/numav.hpp"

#include <cstddef>
#include <tuple>
#include <utility>
#include <memory>
#include <vector>

namespace std
{
    template<>
    struct hash<std::pair<uint64_t, uint64_t>> {
        uint64_t operator()(const std::pair<uint64_t,uint64_t> key) const;
    };

    template<>
    struct hash<std::tuple<uint64_t, uint64_t>> {
        uint64_t operator()(const std::tuple<uint64_t,uint64_t> key) const;
    };

    template<>
    struct hash<std::vector<uint64_t>> {
        uint64_t operator()(const std::vector<uint64_t> vec) const;
    };
}
