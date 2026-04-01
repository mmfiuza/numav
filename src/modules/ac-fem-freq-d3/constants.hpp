// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

namespace numav {

template<ElementOrder O> constexpr size_t NODES_IN_SFC_ELEM = [] {
    if constexpr (O == ElementOrder::O1) { return 3; }
    if constexpr (O == ElementOrder::O2) { return 6; }
    return 0;
}();

template<ElementOrder O> constexpr size_t EXTRA_NODES_IN_SFC_ELEM = [] {
    return NODES_IN_SFC_ELEM<O> - NODES_IN_SFC_ELEM<ElementOrder::O1>;
}();

template<ElementOrder O> constexpr size_t NODES_IN_VOL_ELEM = [] {
    if constexpr (O == ElementOrder::O1) { return 4;  }
    if constexpr (O == ElementOrder::O2) { return 10; }
    return 0;
}();

template<ElementOrder O> constexpr size_t EXTRA_NODES_IN_VOL_ELEM = [] {
    return NODES_IN_VOL_ELEM<O> - NODES_IN_VOL_ELEM<ElementOrder::O1>;
}();

// dimension count in space (not the dimension of mesh elements)
static constexpr size_t DIM = 3; // TODO: review static keyword
static constexpr double PENALTY_METHOD_CONSTANT = 1e12; // TODO: review static

} // namespace numav
