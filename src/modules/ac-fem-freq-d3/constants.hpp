// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "common/aliases.hpp"

namespace numav {

template<ElementOrder O> constexpr size_t NODES_IN_SFC_ELEM = [] {
    if constexpr (O == ElementOrder::O1) { return 3UL; }
    if constexpr (O == ElementOrder::O2) { return 6UL; }
    return 0UL;
}();

template<ElementOrder O> constexpr size_t EXTRA_NODES_IN_SFC_ELEM = [] {
    return NODES_IN_SFC_ELEM<O> - NODES_IN_SFC_ELEM<ElementOrder::O1>;
}();

template<ElementOrder O> constexpr size_t NODES_IN_VOL_ELEM = [] {
    if constexpr (O == ElementOrder::O1) { return 4UL;  }
    if constexpr (O == ElementOrder::O2) { return 10UL; }
    return 0UL;
}();

template<ElementOrder O> constexpr size_t EXTRA_NODES_IN_VOL_ELEM = [] {
    return NODES_IN_VOL_ELEM<O> - NODES_IN_VOL_ELEM<ElementOrder::O1>;
}();

// dimension count in space (not the dimension of mesh elements)
constexpr size_t DIM = 3UL;
constexpr _cmplx_t PENALTY_METHOD_CONSTANT = _cmplx_t(1e12, 0.0);
constexpr size_t NUMAV_DEFAULT_FREQ_STEPS_COUNT = 4096UL;

} // namespace numav
