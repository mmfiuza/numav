// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include <stdexcept>

#include "spdlog/spdlog.h"
#include "fmt/core.h"

namespace numav {

    class Logger {
    public:
        Logger() {
            spdlog::set_level(spdlog::level::level_enum::warn);
        }

        template<typename... T>
        void warn(fmt::format_string<T...> fmt, T&&... args) const {
            spdlog::warn(fmt, args...);
        }
        
        template<typename... T>
        void error(fmt::format_string<T...> fmt, T&&... args) const {
            spdlog::error(fmt, args...);
            throw std::runtime_error("numav error");
        }
    };

} // namespace numav
