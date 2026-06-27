// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/impl.hpp"

namespace numav {

template <ElementOrder O>
void SimulationAcFemFreqD3Tet<O>::Impl::_organize_volume_physical_group_data()
{
    _vei_to_ivpg = fz::SafePtr<uint64_t>(_vei_count);
    for (uint64_t vei = 0; vei != _vei_count; ++vei) {
        _vei_to_ivpg[vei] = _evpg_ivpg_bimap.left.at(_vei_to_evpg[vei]);
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3Tet<O>::Impl::_organize_impedance_physical_group_data()
{
    _isei_count = 0UL;
    for (uint64_t sei = 0UL; sei != _sei_count; ++sei) {
        if (_espg_ispgi_bimap.left.count(_sei_to_espg[sei]) > 0) {
            ++_isei_count;
        }
    }
    _isei_to_sei = fz::SafePtr<uint64_t>(_isei_count);
    uint64_t isei = 0UL;
    for (uint64_t sei = 0UL; sei != _sei_count; ++sei) {
        if (_espg_ispgi_bimap.left.count(_sei_to_espg[sei]) > 0) {
            _isei_to_sei[isei] = sei;
            ++isei;
        }
    }
    _isei_to_ispgi = fz::SafePtr<uint64_t>(_isei_count);
    for (uint64_t isei = 0UL; isei != _isei_count; ++isei) {
        const uint64_t sei = _isei_to_sei[isei];
        _isei_to_ispgi[isei] = _espg_ispgi_bimap.left.at(_sei_to_espg[sei]);
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3Tet<O>::Impl::_organize_velocity_physical_group_data()
{
    _vsei_count = 0UL;
    for (uint64_t sei = 0UL; sei != _sei_count; ++sei) {
        if (_espg_ispgv_bimap.left.count(_sei_to_espg[sei]) > 0) {
            ++_vsei_count;
        }
    }
    _vsei_to_sei = fz::SafePtr<uint64_t>(_vsei_count);
    uint64_t vsei = 0UL;
    for (uint64_t sei = 0UL; sei != _sei_count; ++sei) {
        if (_espg_ispgv_bimap.left.count(_sei_to_espg[sei]) > 0) {
            _vsei_to_sei[vsei] = sei;
            ++vsei;
        }
    }
    _vsei_to_ispgv = fz::SafePtr<uint64_t>(_vsei_count);
    for (uint64_t vsei = 0UL; vsei != _vsei_count; ++vsei) {
        const uint64_t sei = _vsei_to_sei[vsei];
        _vsei_to_ispgv[vsei] = _espg_ispgv_bimap.left.at(_sei_to_espg[sei]);
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3Tet<O>::Impl::_organize_pressure_physical_group_data()
{
    _psei_count = 0UL;
    for (uint64_t sei = 0UL; sei != _sei_count; ++sei) {
        if (_espg_ispgp_bimap.left.count(_sei_to_espg[sei]) > 0) {
            ++_psei_count;
        }
    }
    _psei_to_sei = fz::SafePtr<uint64_t>(_psei_count);
    uint64_t psei = 0UL;
    for (uint64_t sei = 0UL; sei != _sei_count; ++sei) {
        if (_espg_ispgp_bimap.left.count(_sei_to_espg[sei]) > 0) {
            _psei_to_sei[psei] = sei;
            ++psei;
        }
    }
    _psei_to_ispgp = fz::SafePtr<uint64_t>(_psei_count);
    for (uint64_t psei = 0UL; psei != _psei_count; ++psei) {
        const uint64_t sei = _psei_to_sei[psei];
        _psei_to_ispgp[psei] = _espg_ispgp_bimap.left.at(_sei_to_espg[sei]);
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3Tet<O>::Impl::_organize_physical_group_data() {
    _organize_volume_physical_group_data();
    _organize_impedance_physical_group_data();
    _organize_velocity_physical_group_data();
    _organize_pressure_physical_group_data();
}

} // namespace numav

NUMAV_INSTANTIATE_SIM_AC_FEM_FREQ_D3
