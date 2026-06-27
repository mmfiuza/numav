// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include <memory>
#include <array>
#include <vector>

#include "numav/aliases.hpp"

namespace numav {

    enum class Phenomenon : size_t {
        ACOUSTIC
    };

    enum class NumericalMethod : size_t {
        FEM
    };

    enum class Domain : size_t {
        FREQUENCY
    };

    enum class Dimension : size_t {
        D3
    };

    enum class SourceType : size_t {
        POINT,
        SURFACE
    };

    enum class PhysicalQuantity : size_t {
        VOLUME_VELOCITY,
        PARTICLE_VELOCITY,
        PRESSURE,
        IMPEDANCE
    };

    // declare the general Simulation class
    template<
        Phenomenon PHE,
        NumericalMethod NUM,
        Domain DOM,
        Dimension DIM,
        auto... EXTRA
    >
    class Simulation {};

} // namespace numav

// include modules
#include "numav/modules/ac-fem-freq-d3.hpp"
