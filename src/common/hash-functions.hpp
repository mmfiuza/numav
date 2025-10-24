// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "common/debug-macros.hpp"
#include "common/aliases.hpp"

#include <cstddef>
#include <tuple>
#include <utility>
#include <memory>

namespace std
{
    template<>
    struct hash<std::pair<size_t, size_t>> {
        size_t operator()(const std::pair<size_t,size_t>& key) const;
    };

    template<>
    struct hash<std::tuple<size_t, size_t>> {
        size_t operator()(const std::tuple<size_t,size_t>& key) const;
    };

    template<>
    struct hash<std::vector<size_t>> {
        size_t operator()(const std::vector<size_t>& vec) const;
    };
}
