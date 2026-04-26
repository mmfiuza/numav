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

    enum class ElementOrder {
        O1,
        O2
    };

    enum class TypeOfSource {
        POINT,
        SURFACE
    };

    enum class PhysicalQuantity {
        PRESSURE,
        VOLUME_VELOCITY,
        PARTICLE_VELOCITY
    };

    // declare the general Simulation class
    template<
        Phenomenon PHENOMENON,
        NumericalMethod NUMERICAL_METHOD,
        Domain DOMAIN,
        Dimension DIMENSION,
        auto... EXTRA
    >
    class Simulation {};

} // namespace numav

// include modules
#include "numav/modules/ac-fem-freq-d3.hpp"
