// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav/numav.hpp"
#include "common/log.hpp"

namespace numav {

template<ElementOrder O>
SimulationFemHelmTet<O>::Simulation() {
    log::set_level();
    log::set_pattern();
    _is_freq_defined = false;
    _is_mesh_defined = false;
    _is_any_source_defined = false;
    _did_run = false;
    _fi_count = DEFAULT_FREQ_STEPS_COUNT;
    _ivpg_count = 0UL;
    _ispgi_count = 0UL;
    _vpi_count = 0UL;
    _ispgv_count = 0UL;
    _ppi_count = 0UL;
    _ispgp_count = 0UL;
    #if NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
        _solver = std::make_unique<Eigen::SparseLU<
            Eigen::SparseMatrix<Cmplx, Eigen::ColMajor, Eigen::Index>,
            Eigen::COLAMDOrdering<Eigen::Index>
        >>();
    #endif
}

template<ElementOrder O>
SimulationFemHelmTet<O>::~Simulation() {
    if (!_did_run && _is_mesh_defined) {
        _ni_to_coords.free();
        _sei_to_ni.free();
        _vei_to_ni.free();
        _sei_to_espg.free();
        _vei_to_evpg.free();
    }
    if (!_did_run && _is_freq_defined) {
        _fi_to_freq.free();
    }
    if (_did_run) {
        _fi_to_freq.free();
        _a_vals.free();
        _b_row_idx.free();
        _b_vals.free();
        for (uint64_t ivpg = 0UL; ivpg != _ivpg_count; ++ivpg) {
            _ivpg_to_stif_fi_part[ivpg].free();
            _ivpg_to_mass_fi_part[ivpg].free();
            _ivpg_to_ptr_in_a[ivpg].free();
        }
        _ivpg_to_stif_fi_part.free();
        _ivpg_to_mass_fi_part.free();
        _ivpg_to_ptr_in_a.free();
        for (uint64_t ispgi = 0UL; ispgi != _ispgi_count; ++ispgi) {
            _ispgi_to_damp_fi_part[ispgi].free();
            _ispgi_to_ptr_in_a[ispgi].free();
        }
        _ispgi_to_damp_fi_part.free();
        _ispgi_to_ptr_in_a.free();
        _vpi_to_ptr_in_b.free();
        for (uint64_t ispgv = 0UL; ispgv != _ispgv_count; ++ispgv) {
            _ispgv_to_forc_fi_part[ispgv].free();
            _ispgv_to_ptr_in_b[ispgv].free();
        }
        _ispgv_to_forc_fi_part.free();
        _ispgv_to_ptr_in_b.free();
        _apvi_to_pressure.free();
        for (uint64_t apvi = 0UL; apvi != _apvi_count; ++apvi) {
            _apvi_to_ptr_in_a[apvi].free();
            _apvi_to_ptr_in_b[apvi].free();
        }
        _apvi_to_ptr_in_a.free();
        _apvi_to_ptr_in_b.free();
        _x.free();
        #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
            _b_dense.free();
            _a_diag.free();
        #elif NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
            _a_row_idx.free();
            _a_col_idx.free();
            _b_row_idx_signed.free();
        #elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
            _b_dense.free();
        #elif NUMAV_SYSTEM_SOLVER == NUMAV_MUMPS
            _b_dense.free();
        #endif
    }
}

template<ElementOrder O>
SimulationFemHelmTet<O>::Simulation(
    Simulation&& other
) noexcept = default;

template<ElementOrder O>
SimulationFemHelmTet<O>& SimulationFemHelmTet<O>::operator=(
    Simulation&& other
) noexcept = default;

} // namespace numav

NUMAV_INSTANTIATE_SIM_AC_FEM_FREQ_D3
