// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/impl.hpp"

#include "common/log.hpp"

namespace numav {

template <ElementOrder O>
SimulationAcFemFreqD3<O>::Simulation() {
    pimpl = std::make_unique<Impl>();
}

template <ElementOrder O>
SimulationAcFemFreqD3<O>::~Simulation() = default;

template<ElementOrder O>
SimulationAcFemFreqD3<O>::Impl::Impl() {
    log::set_level();
    log::set_pattern();
    _is_mesh_defined = false;
    _is_freq_range_defined = false;
    _is_any_source_defined = false;
    _did_run = false;
}

template<ElementOrder O>
SimulationAcFemFreqD3<O>::Impl::~Impl() {
    _freq_steps.free();
    _node_coords.free();
    _sei_to_ni.free();
    _vei_to_ni.free();
    _sei_to_espg.free();
    _vei_to_evpg.free();
    _nnz_rowcol_idx_pairs.free();
    for (size_t ispgv=0; ispgv!=_ispgv_count(); ++ispgv) {
        _ispgv_to_forc_fi_part[ispgv].free();
        _ispgv_to_ptr_in_b[ispgv].free();
    }
    _pvni_to_forc_fi_part.free();
    _pvni_to_ptr_in_b.free();
    _ispgv_to_forc_fi_part.free();
    _ispgv_to_ptr_in_b.free();
    for (size_t ivpg=0; ivpg!=_ivpg_count(); ++ivpg) {
        _ivpg_to_stif_fi_part[ivpg].free();
        _ivpg_to_mass_fi_part[ivpg].free();
        _ivpg_to_ptr_in_a[ivpg].free();
    }
    _ivpg_to_stif_fi_part.free();
    _ivpg_to_mass_fi_part.free();
    _ivpg_to_ptr_in_a.free();
    for (size_t ispgi=0; ispgi!=_ispgi_count(); ++ispgi) {
        _ispgi_to_damp_fi_part[ispgi].free();
        _ispgi_to_ptr_in_a[ispgi].free();
    }
    _ispgi_to_damp_fi_part.free();
    _ispgi_to_ptr_in_a.free();
    _isei_to_ispgi.free();
    _vsei_to_ispgv.free();
    _vei_to_ivpg.free();
    _ispgv_to_velocity.free();
    _ispgp_to_pressure.free();
    _ispgi_to_impedance.free();
    _ivpg_to_volprop.free();
    _a_vals.free();
    _b_row_idx.free();
    _b_vals.free();
    _isei_to_sei.free();
    _vsei_to_sei.free();
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
