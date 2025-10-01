// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "spdlog/spdlog.h"
#include "fmt/core.h"

#include <stdexcept>

namespace numav::log {

    void set_level();

    void set_pattern();

    template<typename... T>
    void warn(fmt::format_string<T...> fmt, T&&... args) {
        spdlog::warn(fmt, args...);
    }

    template<typename... T>
    void error(fmt::format_string<T...> fmt, T&&... args) {
        spdlog::error(fmt, args...);
        throw std::runtime_error("numav error");
    }

} // namespace numav::log
