// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include <stdexcept>

#include "spdlog/spdlog.h"
#include "fmt/core.h"

namespace numav::log {

    void set_level() {
        spdlog::set_level(spdlog::level::level_enum::warn);
    }

    void set_pattern() {
        spdlog::set_pattern("[%^numav %l%$]: %v");
    }

    template<typename... T>
    void warn(fmt::format_string<T...> fmt, T&&... args) {
        spdlog::warn(fmt, args...);
    }

    template<typename... T>
    void error(fmt::format_string<T...> fmt, T&&... args) {
        spdlog::error(fmt, args...);
        throw std::runtime_error("numav error");
    }

} // namespace numav
