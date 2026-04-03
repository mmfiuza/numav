// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/impl.hpp"

#include "common/maths.hpp"

namespace numav {

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_define_freq_vector() {
    // todo: decide number here
    // todo: make it not linear
    _freq_count = 4096;
    _freq_steps = cubespace(_freq_min, _freq_max, _freq_count);
}

// explicit instantiation declarations
INSTANTIATE_SIMULATION_CLASS

} // namespace numav