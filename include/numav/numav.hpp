// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include <memory>
#include <array>
#include <vector>
#include <functional>
#include <complex>

namespace numav {

    // type aliases
    using Float = double;
    inline constexpr Float operator"" _F(unsigned long long v) noexcept {
        return static_cast<Float>(v);
    }
    inline constexpr Float operator"" _F(long double v) noexcept {
        return static_cast<Float>(v);
    }
    using Cmplx = typename std::complex<Float>;
    using FuncFloatToCmplx = typename std::function<Cmplx(const Float)>;
    using Coord = typename std::array<Float,3UL>;

    // constants
    constexpr Float PI = 3.14159265358979323846_F;

    // numerical methods
    enum class NumericalMethod : uint64_t {
        FEM
    };

    // declare the general Simulation class
    template<NumericalMethod N, auto... ARGS> class Simulation {};

} // namespace numav

// include modules
#include "numav/modules/fem/fem.hpp"
