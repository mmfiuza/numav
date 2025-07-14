// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include <array>
#include <complex>
#include <functional>
#include <unordered_map>

#define SAFE_PTR_DEBUG
#include "SafePtr.hpp"

#include "Eigen/Dense"

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

    template<Dimension D> constexpr size_t DIM_COUNT = []{
        if constexpr (D == Dimension::D1) return 1;
        if constexpr (D == Dimension::D2) return 2;
        if constexpr (D == Dimension::D3) return 3;
        return 0;
    }();

    enum class ElementOrder {
        O1,
        O2
    };

    template<ElementOrder O> constexpr size_t NODES_IN_2D_ELEM = []{
        if constexpr (O == ElementOrder::O1) return 3;
        if constexpr (O == ElementOrder::O2) return 6;
        return 0;
    }();

    template<ElementOrder O> constexpr size_t EXTRA_NODES_IN_2D_ELEM = []{
        return NODES_IN_2D_ELEM<O> - NODES_IN_2D_ELEM<ElementOrder::O1>;
    }();

    template<ElementOrder O> constexpr size_t NODES_IN_3D_ELEM = []{
        if constexpr (O == ElementOrder::O1) return 4;
        if constexpr (O == ElementOrder::O2) return 10;
        return 0;
    }();

    template<ElementOrder O> constexpr size_t EXTRA_NODES_IN_3D_ELEM = []{
        return NODES_IN_3D_ELEM<O> - NODES_IN_3D_ELEM<ElementOrder::O1>;
    }();

    enum class TypeOfSource {
        POINT,
        SURFACE,
    };

    enum class PhysicalQuantity {
        PRESSURE,
        VOLUME_VELOCITY,
    };

    // declare the general Result class
    template<
        Phenomenon PHENOMENON,
        NumericalMethod NUMERICAL_METHOD,
        Domain DOMAIN,
        Dimension DIMENSION
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

// include specific numerical methods
#include "ac-fem-freq-d3.hpp"
