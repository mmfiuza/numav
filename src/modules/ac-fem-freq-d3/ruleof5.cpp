// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/impl.hpp"

#include "common/log.hpp"

namespace numav {

template <ElementOrder O>
SimulationAcFemFreqD3<O>::Simulation() {
    _pimpl = std::make_unique<Impl>();
}

template <ElementOrder O>
SimulationAcFemFreqD3<O>::~Simulation() = default;

template<ElementOrder O>
SimulationAcFemFreqD3<O>::Simulation(
    Simulation&& other
) noexcept = default;

template<ElementOrder O>
SimulationAcFemFreqD3<O>& SimulationAcFemFreqD3<O>::operator=(
    Simulation&& other
) noexcept = default;

template<ElementOrder O>
SimulationAcFemFreqD3<O>::Impl::Impl() {
    log::set_level();
    log::set_pattern();
    _is_mesh_defined = false;
    _is_any_source_defined = false;
    _did_run = false;
    _freq_type_defined_by_user = _FreqTypeDefinedByUser::UNDEFINED;
    _frequency_sampling_density = FrequencySamplingDensity::QUADRATIC;
    _freq_count = NUMAV_DEFAULT_FREQ_STEPS_COUNT;
    _ivpg_count = 0;
    _ispgi_count = 0;
    _pvni_count = 0;
    _ispgv_count = 0;
    _ppni_count = 0;
    _ispgp_count = 0;
    _ri_count = 0;
}

template<ElementOrder O>
SimulationAcFemFreqD3<O>::Impl::~Impl() {
    _freq_steps.free();
    _ni_to_coords.free();
    _sei_to_ni.free();
    _vei_to_ni.free();
    _sei_to_espg.free();
    _vei_to_evpg.free();
    _isei_to_sei.free();
    _vsei_to_sei.free();
    _psei_to_sei.free();
    _vei_to_ivpg.free();
    _isei_to_ispgi.free();
    _vsei_to_ispgv.free();
    _psei_to_ispgp.free();
    _nnz_rowcol_idx_pairs.free();
    _a_vals.free();
    _b_row_idx.free();
    _b_vals.free();
    _ivpg_to_volprop.free();
    _ispgi_to_impedance.free();
    _ispgv_to_velocity.free();
    _ispgp_to_pressure.free();
    for (size_t ivpg=0; ivpg!=_ivpg_count; ++ivpg) {
        _ivpg_to_stif_fi_part[ivpg].free();
        _ivpg_to_mass_fi_part[ivpg].free();
        _ivpg_to_ptr_in_a[ivpg].free();
    }
    _ivpg_to_stif_fi_part.free();
    _ivpg_to_mass_fi_part.free();
    _ivpg_to_ptr_in_a.free();
    for (size_t ispgi=0; ispgi!=_ispgi_count; ++ispgi) {
        _ispgi_to_damp_fi_part[ispgi].free();
        _ispgi_to_ptr_in_a[ispgi].free();
    }
    _ispgi_to_damp_fi_part.free();
    _ispgi_to_ptr_in_a.free();
    _pvni_to_forc_fi_part.free();
    _pvni_to_ptr_in_b.free();
    for (size_t ispgv=0; ispgv!=_ispgv_count; ++ispgv) {
        _ispgv_to_forc_fi_part[ispgv].free();
        _ispgv_to_ptr_in_b[ispgv].free();
    }
    _ispgv_to_forc_fi_part.free();
    _ispgv_to_ptr_in_b.free();
    _pvi_to_pressure.free();
    for (size_t pvi=0; pvi!=_pvi_count; ++pvi) {
        _pvi_to_ptr_in_a[pvi].free();
        _pvi_to_ptr_in_b[pvi].free();
    }
    _pvi_to_ptr_in_a.free();
    _pvi_to_ptr_in_b.free();
    _x.free();
    #if NUMAV_SYSTEM_SOLVER == NUMAV_LDLT_SOLVER
        _b_dense.free();
        _a_diag.free();
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
        _b_dense.free();
    #endif
}

template<ElementOrder O>
SimulationAcFemFreqD3<O>::Impl::Impl(
    Impl&& other
) noexcept = default;

template<ElementOrder O>
SimulationAcFemFreqD3<O>::Impl& SimulationAcFemFreqD3<O>::Impl::operator=(
    Impl&& other
) noexcept = default;

// explicit instantiation declarations
INSTANTIATE_SIMULATION_CLASS

} // namespace numav
