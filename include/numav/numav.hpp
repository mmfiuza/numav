// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

/*
Glossary:

    entity:
        Some geometrical body (points, lines, surfaces, volumes, etc).

    index (idx):
        An unsigned integer that uniquely identifies an entity.

    EPG (External Physical Group):
        An unsigned integer present in an external mesh file (.bdf, etc) or
        API that classifies an entity into some group, usually a material. An
        EPG sequence is NOT necessarily consecutive, since mesh formats can
        have arbitrary conventions.

    IPG (Internal Physical Group)
        An unsigned integer that classifies an entity into some group in Numav's
        internal implementation. An IPG sequence is ALWAYS consecutive starting
        at 0, because that is conveninent to manipulate in C++.
*/

#pragma once

#include <iostream>
#include <array>
#include <complex>
#include <functional>

namespace numav {

    enum class Phenomenon {
        ACOUSTIC,
    };

    enum class NumericalMethod {
        FEM,
    };

    enum class Domain {
        FREQUENCY,
        TIME
    };

    enum class Dimension {
        D1,
        D2,
        D3,
    };

    enum class ElementOrder {
        O1,
        O2
    };

    enum class TypeOfSource {
        POINT,
        SURFACE,
    };

    enum class PhysicalQuantity {
        PRESSURE,
        VOLUME_VELOCITY,
    };
    
    template<Dimension D> constexpr size_t DIM_COUNT = [] {
        if constexpr (D == Dimension::D1) { return 1; }
        if constexpr (D == Dimension::D2) { return 2; }
        if constexpr (D == Dimension::D3) { return 3; }
        return 0;
    }();

    // declare the general Result class
    template<
        Phenomenon PHENOMENON,
        NumericalMethod NUMERICAL_METHOD,
        Domain DOMAIN,
        Dimension DIMENSION,
        auto... EXTRA
    >
    class Result {};

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
