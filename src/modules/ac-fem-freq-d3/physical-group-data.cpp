// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/impl.hpp"

namespace numav {

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_organize_volume_physical_group_data()
{
    _ivpg_to_volprop = fz::SafePtr<_VolProp>(_ivpg_count);
    for (const auto& [evpg, volprop] : _evpg_to_volprop) {
        const size_t ivpg = _evpg_to_ivpg.at(evpg);
        _ivpg_to_volprop[ivpg] = volprop;
    }
    _vei_to_ivpg = fz::SafePtr<size_t>(_vei_count);
    for (size_t vei = 0; vei != _vei_count; ++vei) {
        _vei_to_ivpg[vei] = _evpg_to_ivpg.at(_vei_to_evpg[vei]);
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_organize_impedance_physical_group_data()
{
    _ispgi_to_impedance = fz::SafePtr<FuncFloatToCmplx>(_ispgi_count);
    for (const auto& [espg, impedance] : _espg_to_impedance) {
        const size_t ispgi = _espg_to_ispg.at(espg);
        _ispgi_to_impedance[ispgi] = impedance;
    }
    _isei_count = 0UL;
    for (size_t sei = 0UL; sei != _sei_count; ++sei) {
        if (_espg_to_impedance.contains(_sei_to_espg[sei])) {
            ++_isei_count;
        }
    }
    _isei_to_sei = fz::SafePtr<size_t>(_isei_count);
    size_t isei = 0UL;
    for (size_t sei = 0UL; sei != _sei_count; ++sei) {
        if (_espg_to_impedance.contains(_sei_to_espg[sei])) {
            _isei_to_sei[isei] = sei;
            ++isei;
        }
    }
    _isei_to_ispgi = fz::SafePtr<size_t>(_isei_count);
    for (size_t isei = 0UL; isei != _isei_count; ++isei) {
        const size_t sei = _isei_to_sei[isei];
        _isei_to_ispgi[isei] = _espg_to_ispg.at(_sei_to_espg[sei]);
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_organize_velocity_physical_group_data()
{
    _ispgv_to_velocity = fz::SafePtr<FuncFloatToCmplx>(_ispgv_count);
    for (const auto& [espg, volvel] : _espg_to_velocity) {
        const size_t ispgv = _espg_to_ispg.at(espg);
        _ispgv_to_velocity[ispgv] = volvel;
    }
    _vsei_count = 0UL;
    for (size_t sei = 0UL; sei != _sei_count; ++sei) {
        if (_espg_to_velocity.contains(_sei_to_espg[sei])) {
            ++_vsei_count;
        }
    }
    _vsei_to_sei = fz::SafePtr<size_t>(_vsei_count);
    size_t vsei = 0UL;
    for (size_t sei = 0UL; sei != _sei_count; ++sei) {
        if (_espg_to_velocity.contains(_sei_to_espg[sei])) {
            _vsei_to_sei[vsei] = sei;
            ++vsei;
        }
    }
    _vsei_to_ispgv = fz::SafePtr<size_t>(_vsei_count);
    for (size_t vsei = 0UL; vsei != _vsei_count; ++vsei) {
        const size_t sei = _vsei_to_sei[vsei];
        _vsei_to_ispgv[vsei] = _espg_to_ispg.at(_sei_to_espg[sei]);
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_organize_pressure_physical_group_data()
{
    _ispgp_to_pressure = fz::SafePtr<FuncFloatToCmplx>(_ispgp_count);
    for (const auto& [espg, pressure] : _espg_to_pressure) {
        const size_t ispgp = _espg_to_ispg.at(espg);
        _ispgp_to_pressure[ispgp] = pressure;
    }
    _psei_count = 0UL;
    for (size_t sei = 0UL; sei != _sei_count; ++sei) {
        if (_espg_to_pressure.contains(_sei_to_espg[sei])) {
            ++_psei_count;
        }
    }
    _psei_to_sei = fz::SafePtr<size_t>(_psei_count);
    size_t psei = 0UL;
    for (size_t sei = 0UL; sei != _sei_count; ++sei) {
        if (_espg_to_pressure.contains(_sei_to_espg[sei])) {
            _psei_to_sei[psei] = sei;
            ++psei;
        }
    }
    _psei_to_ispgp = fz::SafePtr<size_t>(_psei_count);
    for (size_t psei = 0UL; psei != _psei_count; ++psei) {
        const size_t sei = _psei_to_sei[psei];
        _psei_to_ispgp[psei] = _espg_to_ispg.at(_sei_to_espg[sei]);
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_organize_physical_group_data() {
    _organize_volume_physical_group_data();
    _organize_impedance_physical_group_data();
    _organize_velocity_physical_group_data();
    _organize_pressure_physical_group_data();
}

} // namespace numav

NUMAV_INSTANTIATE_SIM_AC_FEM_FREQ_D3
