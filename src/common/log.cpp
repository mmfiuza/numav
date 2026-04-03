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
        std::cout << std::string(80,'-') << "\n";
        std::cout << std::string(34,' ') <<
            "Numav Solver" << std::string(34,' ') << "\n";
        std::cout << std::string(20,' ') <<
            "Copyright (c) 2026 Matheus Machado Fiuza" <<
            std::string(20,' ') << "\n";
        std::cout << std::string(30,' ') << "All rights reserved." << 
            std::string(30,' ') << "\n";
        std::cout << std::string(80, '-') << "\n";
    }

    void print_opening_ac_fem_freq_d3() {
        std::cout << "Acoustic 3D frequency domain FEM simulation started.\n";
    }

} // namespace numav::log
