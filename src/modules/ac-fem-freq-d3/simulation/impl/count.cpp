// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav/numav.hpp"
#include "modules/ac-fem-freq-d3/simulation/impl/impl.hpp"

namespace numav {

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::Impl::_ni_count() const {
    return _node_coords.size();
}

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::Impl::_sei_count() const {
    return _sei_to_ni.size();
}

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::Impl::_vei_count() const {
    return _vei_to_ni.size();
}

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::Impl::_isei_count() const {
    return _isei_to_ispgi.size();
}

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::Impl::_vsei_count() const {
    return _vsei_to_ispgv.size();
}

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::Impl::_pvni_count() const {
    return _point_volvel.size();
}

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::Impl::_ispgi_count() const {
    return _espg_to_impedance.size();
}

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::Impl::_ispgv_count() const {
    return _espg_to_velocity.size();
}

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::Impl::_ivpg_count() const {
    return _evpg_to_volprop.size();
}

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::Impl::_freq_count() const {
    return _freq_steps.size();
}

// explicit instantiation declarations
template class Simulation<
    Phenomenon::ACOUSTIC,
    NumericalMethod::FEM,
    Domain::FREQUENCY,
    Dimension::D3,
    ElementOrder::O1
>;
template class Simulation<
    Phenomenon::ACOUSTIC,
    NumericalMethod::FEM,
    Domain::FREQUENCY,
    Dimension::D3,
    ElementOrder::O2
>;

} // namespace numav
