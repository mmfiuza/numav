// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/log.hpp"

#include <iostream>
#include <iomanip>

#include <indicators/cursor_control.hpp>

namespace numav::log {

    void set_level() {
        spdlog::set_level(spdlog::level::level_enum::warn);
    }

    void set_pattern() {
        spdlog::set_pattern("[%^numav %l%$]: %v");
    }

    void print_opening() {
        std::cout << std::string(80UL,'-') << "\n";
        std::cout << std::string(37UL,' ') <<
            "Numav" << std::string(38UL,' ') << "\n";
        std::cout << std::string(20UL,' ') <<
            "Copyright (c) 2026 Matheus Machado Fiuza" <<
            std::string(20UL,' ') << "\n";
        std::cout << std::string(21UL,' ') << 
            "Licensed under the AGPL3 Licence." << 
            std::string(22UL,' ') << "\n";
        std::cout << std::string(80UL, '-') << "\n";
    }

    void print_opening_ac_fem_freq_d3() {
        std::cout << "Acoustic 3D frequency domain FEM simulation started.\n";
    }

    void print_start_time() {
        const auto time = std::chrono::system_clock::now();
        const time_t time_time_t = std::chrono::system_clock::to_time_t(time);
        std::cout <<
            "Solver started at: " <<
            std::put_time(std::localtime(&time_time_t), "%Hh:%Mm:%Ss") << 
            "\n";
    }

    void print_finish_time() {
        const auto time = std::chrono::system_clock::now();
        const time_t time_time_t = std::chrono::system_clock::to_time_t(time);
        std::cout <<
            "Solver finished at: " <<
            std::put_time(std::localtime(&time_time_t), "%Hh:%Mm:%Ss") <<
            "\n";
    }

    size_t start_progress_bar(
        std::unique_ptr<indicators::ProgressBar>& bar,
        const size_t progress_max
    ) {
        indicators::show_console_cursor(false);
        bar = std::make_unique<indicators::ProgressBar>(
            indicators::option::BarWidth{37UL},
            indicators::option::Start{" |"},
            indicators::option::Fill{"="},
            indicators::option::Lead{"="},
            indicators::option::Remainder{" "},
            indicators::option::End{"|"},
            indicators::option::PrefixText{"Running"},
            indicators::option::ForegroundColor{indicators::Color::unspecified},
            indicators::option::ShowPercentage{true},
            indicators::option::ShowElapsedTime{true},
            indicators::option::ShowRemainingTime{true},
            indicators::option::FontStyles{
                std::vector<indicators::FontStyle>{indicators::FontStyle::bold}
            },
            indicators::option::MinProgress{0UL},
            indicators::option::MaxProgress{progress_max}
        );
        size_t progress = 0UL;
        bar->set_progress(progress);
        return progress;
    }

    void increment_progress_bar(
        std::unique_ptr<indicators::ProgressBar>& bar,
        size_t& bar_progress
    ) {
        ++bar_progress;
        bar->set_progress(bar_progress);
    }

    void finish_progress_bar() {
        indicators::show_console_cursor(true);
    }

} // namespace numav::log
