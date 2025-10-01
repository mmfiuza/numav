// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/log.hpp"

namespace numav::log {

    void set_level() {
        spdlog::set_level(spdlog::level::level_enum::warn);
    }

    void set_pattern() {
        spdlog::set_pattern("[%^numav %l%$]: %v");
    }

} // namespace numav::log
