// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/debug-macros.hpp"
#include "numav/aliases.hpp"

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
}
