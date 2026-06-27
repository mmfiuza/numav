// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include <memory>
#include <array>
#include <vector>

#include "numav/aliases.hpp"

namespace numav {

    enum class Phenomenon : uint64_t {
        ACOUSTIC
    };

    enum class NumericalMethod : uint64_t {
        FEM
    };

    enum class Domain : uint64_t {
        FREQUENCY
    };

    enum class Dimension : uint64_t {
        D3
    };

    enum class SourceType : uint64_t {
        POINT,
        SURFACE
    };

    enum class PhysicalQuantity : uint64_t {
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
