// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/aliases.hpp"
#include "common/log.hpp"

#include "fmt/core.h"

#include <stdexcept>

namespace numav {

    template<typename... T>
    void error(fmt::format_string<T...> fmt, T&&... args) {
        log::error(fmt, args...);
        throw std::runtime_error("numav error");
    }

} // namespace numav
