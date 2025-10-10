// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

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
        std::cout << "----------------------------------------\n";
        std::cout << "Numav Solver\n";
        std::cout << "Copyright (c) 2025 Numav\n";
        std::cout << "All rights reserved.\n";
        std::cout << "----------------------------------------\n";
    }

    void print_opening_ac_fem_freq_d3() {
        std::cout << "Acoustic 3D frequency domain FEM simulation started.\n";
    }

} // namespace numav::log
