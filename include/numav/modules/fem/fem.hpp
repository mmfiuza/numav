// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

namespace numav {

constexpr uint64_t DIM = 3UL; // dimension count in space
constexpr Cmplx PENALTY_METHOD_CONSTANT = Cmplx(1e12_F, 0_F);

enum class Equation : uint64_t {
    HELMHOLTZ
};

enum class ElementShape : uint64_t {
    TETRAHEDRON
};

enum class ElementOrder : uint64_t {
    LINEAR,
    QUADRATIC
};

template<ElementOrder O> constexpr uint64_t ENIS_COUNT = [] {
    if constexpr (O == ElementOrder::LINEAR) { return 3UL; }
    if constexpr (O == ElementOrder::QUADRATIC) { return 6UL; }
    return 0UL;
}();

template<ElementOrder O> constexpr uint64_t EXTRA_ENIS_COUNT = [] {
    return ENIS_COUNT<O> - ENIS_COUNT<ElementOrder::LINEAR>;
}();

template<ElementOrder O> constexpr uint64_t ENIV_COUNT = [] {
    if constexpr (O == ElementOrder::LINEAR) { return 4UL;  }
    if constexpr (O == ElementOrder::QUADRATIC) { return 10UL; }
    return 0UL;
}();

template<ElementOrder O> constexpr uint64_t EXTRA_ENIV_COUNT = [] {
    return ENIV_COUNT<O> - ENIV_COUNT<ElementOrder::LINEAR>;
}();

} // namespace numav

// include modules
#include "numav/modules/fem/helmholtz/helmholtz.hpp"
