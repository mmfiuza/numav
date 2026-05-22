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
        FREQUENCY,
        TIME
    };

    enum class Dimension : size_t {
        D1,
        D2,
        D3
    };

    enum class TypeOfSource : size_t { // TODO: rename to SourceType
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
