// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/aliases.hpp"

namespace numav {

    enum class Phenomenon {
        ACOUSTIC
    };

    enum class NumericalMethod {
        FEM
    };

    enum class Domain {
        FREQUENCY,
        TIME
    };

    enum class Dimension {
        D1,
        D2,
        D3
    };

    enum class TypeOfSource {
        POINT,
        SURFACE
    };

    enum class PhysicalQuantity {
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
