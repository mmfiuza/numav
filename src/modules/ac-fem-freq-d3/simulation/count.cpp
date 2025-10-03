// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav/numav.hpp"

namespace numav {

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::_ni_count() const {
    return _node_coords.size();
}

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::_sei_count() const {
    return _sei_to_ni.size();
}

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::_vei_count() const {
    return _vei_to_ni.size();
}

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::_isei_count() const {
    return _isei_to_ispgi.size();
}

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::_ispgi_count() const {
    return _espg_to_impedance.size();
}

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::_ivpg_count() const {
    return _evpg_to_volprop.size();
}

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::_freq_count() const {
    return _freq_steps.size();
}

} // namespace numav
