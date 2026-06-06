// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/aliases.hpp"
#include "common/utils.hpp"

namespace numav {

template<ElementOrder O> constexpr size_t ENI_COUNT_SFC = [] {
    if constexpr (O == ElementOrder::O1) { return 3UL; }
    if constexpr (O == ElementOrder::O2) { return 6UL; }
    return 0UL;
}();

template<ElementOrder O> constexpr size_t EXTRA_ENI_COUNT_SFC = [] {
    return ENI_COUNT_SFC<O> - ENI_COUNT_SFC<ElementOrder::O1>;
}();

template<ElementOrder O> constexpr size_t ENI_COUNT_VOL = [] {
    if constexpr (O == ElementOrder::O1) { return 4UL;  }
    if constexpr (O == ElementOrder::O2) { return 10UL; }
    return 0UL;
}();

template<ElementOrder O> constexpr size_t EXTRA_ENI_COUNT_VOL = [] {
    return ENI_COUNT_VOL<O> - ENI_COUNT_VOL<ElementOrder::O1>;
}();

// dimension count in space (not the dimension of mesh elements)
constexpr size_t DIM = 3UL;
constexpr Cmplx PENALTY_METHOD_CONSTANT = Cmplx(1e12_F, 0_F);
constexpr size_t DEFAULT_FREQ_STEPS_COUNT = 4096UL;
constexpr Float PI = 3.14159265358979323846_F;

} // namespace numav
