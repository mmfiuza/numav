// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/aliases.hpp"

#include "spdlog/spdlog.h"
#include "fmt/core.h"
#include <indicators/progress_bar.hpp>

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
    }

    void print_opening();
    void print_opening_ac_fem_freq_d3();
    void print_start_time();
    void print_finish_time();

    size_t start_progress_bar(
        std::unique_ptr<indicators::ProgressBar>& bar,
        const size_t progress_max
    );
    void increment_progress_bar(
        std::unique_ptr<indicators::ProgressBar>& bar,
        size_t& bar_progress
    );
    void finish_progress_bar();

} // namespace numav::log
