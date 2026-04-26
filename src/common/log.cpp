// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/log.hpp"

#include <iostream>

namespace numav::log {

    void set_level() {
        spdlog::set_level(spdlog::level::level_enum::warn);
    }

    void set_pattern() {
        spdlog::set_pattern("[%^numav %l%$]: %v");
    }

    void print_opening() {
        std::cout << std::string(80UL,'-') << "\n";
        std::cout << std::string(34UL,' ') <<
            "Numav Solver" << std::string(34UL,' ') << "\n";
        std::cout << std::string(20UL,' ') <<
            "Copyright (c) 2026 Matheus Machado Fiuza" <<
            std::string(20UL,' ') << "\n";
        std::cout << std::string(30UL,' ') << "All rights reserved." << 
            std::string(30UL,' ') << "\n";
        std::cout << std::string(80UL, '-') << "\n";
    }

    void print_opening_ac_fem_freq_d3() {
        std::cout << "Acoustic 3D frequency domain FEM simulation started.\n";
    }

} // namespace numav::log
