// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/numav.hpp"
#include "common/log.hpp"

#include <stdexcept>

#include "fmt/core.h"

namespace numav {

    template<typename... T>
    void error(fmt::format_string<T...> fmt, T&&... args) {
        log::error(fmt, args...);
        throw std::runtime_error("numav error");
    }

} // namespace numav
