// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav/numav.hpp"

#include "common/log.hpp"

#include <tuple>
#include <fstream>

namespace numav {

template<ElementOrder O>
SimulationAcFemFreqD3<O>::Simulation() {
    log::set_level();
    log::set_pattern();
    _is_mesh_defined = false;
    _is_freq_range_defined = false;
    _is_any_source_defined = false;
}

template<ElementOrder O>
SimulationAcFemFreqD3<O>::~Simulation() {
    _freq_steps.free();
    _node_coords.free();
    _sei_to_ni.free();
    _vei_to_ni.free();
    _sei_to_espg.free();
    _vei_to_evpg.free();
    _nnz_rowcol_idx_pairs.free();
    for (_ipg_t ivpg=0; ivpg!=_ivpg_count(); ++ivpg) {
        _ivpg_to_stif_fi_part[ivpg].free();
        _ivpg_to_mass_fi_part[ivpg].free();
        _ivpg_to_ptr_in_a[ivpg].free();
    }
    _ivpg_to_stif_fi_part.free();
    _ivpg_to_mass_fi_part.free();
    for (_ipg_t ispgi=0; ispgi!=_ispgi_count(); ++ispgi) {
        _ispgi_to_damp_fi_part[ispgi].free();
        _ispgi_to_ptr_in_a[ispgi].free();
    }
    _ispgi_to_damp_fi_part.free();
    _ispgi_to_ptr_in_a.free();
    _isei_to_ispgi.free();
    _vei_to_ivpg.free();
    _ispgv_to_volvel.free();
    _ispgp_to_pressure.free();
    _ispgi_to_impedance.free();
    _ivpg_to_volprop.free();
    _ivpg_to_ptr_in_a.free();
    _a_vals.free();
    _isei_to_sei.free();
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::_check_if_mesh_is_defined() {
    if (!_is_mesh_defined){
        log::error("Mesh not defined. Call load_mesh to do so.");
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::set_frequency_range(
    const double& freq_min,
    const double& freq_max
) {
    if (freq_min<0 || freq_max<0) {
        log::error("Frequency limits should be positive.");
    }
    if (freq_min >= freq_max) {
        log::error("Upper frequency should be greater than the lower.");
    }
    _freq_min = freq_min;
    _freq_max = freq_max;
    _is_freq_range_defined = true;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::load_mesh(
    const char* const path_to_mesh
) {
    if (_is_mesh_defined) {
        log::error("Mesh is already defined.");
    }
    std::string str = path_to_mesh;
    if (str.ends_with(".bdf") || str.ends_with(".nas")) {
        _load_bdf(path_to_mesh);
    }
    else {
        const size_t dot_position = str.find_last_of('.');
        const size_t format_len = str.size() - dot_position;
        const std::string format = std::string(
            str.substr(dot_position, format_len)
        );
        log::error("Unrecognized file format: \"{0}\".", format);
    }
    _generate_extra_nodes(); // call is based on the element order
    _is_mesh_defined = true;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_volume_material(
    const _epg_t& evpg,
    const _FuncRealToCmplx& density,
    const _FuncRealToCmplx& soundspeed
) {
    _check_if_mesh_is_defined();
    if (!_existing_evpg.contains(evpg)) {
        log::error("Tag {} not found in mesh file.", evpg);
    }
    if (_evpg_to_volprop.contains(evpg)) {
        log::error("Tag {} already assigned.", evpg);
    }
    _evpg_to_volprop.insert({evpg, {density, soundspeed}});
    const _ipg_t ivpg = _evpg_to_ivpg.size();
    _evpg_to_ivpg.insert({evpg, ivpg});
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_sound_source(
    const TypeOfSource& type_of_source,
    const std::array<double,3>& point_coordinates,
    const PhysicalQuantity& physical_quantity_type,
    const _FuncRealToCmplx& physical_quantity_value
) {
    _check_if_mesh_is_defined();
    if (type_of_source != TypeOfSource::POINT) {
        log::error("Tried to assign coordinates to a surface sound source.");
    }
    const _idx_t closest_point_idx = _get_closest_point(point_coordinates);
    
    if (physical_quantity_type == PhysicalQuantity::VOLUME_VELOCITY) {
        _point_volvel.push_back(
            std::make_tuple(closest_point_idx, physical_quantity_value)
        );
    }
    else if (physical_quantity_type == PhysicalQuantity::PRESSURE) {
        _point_pressure.push_back(
            std::make_tuple(closest_point_idx, physical_quantity_value)
        );
    }
    else {
        log::error("Possible physical quantities are volume velocity"
            " or pressure.");
    }
    _is_any_source_defined = true;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_sound_source(
    const TypeOfSource& type_of_source,
    const _epg_t& espg,
    const PhysicalQuantity& physical_quantity_type,
    const _FuncRealToCmplx& physical_quantity_value
) {
    _check_if_mesh_is_defined();
    if (type_of_source != TypeOfSource::SURFACE) {
        log::error("Tried to assign a tag to a point.");
    }
    if (!_existing_espg.contains(espg)) {
        log::error("Tag {} not found in mesh file.", espg);
    }
    if (_espg_to_pressure.contains(espg) ||
        _espg_to_volvel.contains(espg) ||
        _espg_to_impedance.contains(espg)
    ) {
        log::error("Tag {} already assigned.", espg);
    }
    if (physical_quantity_type == PhysicalQuantity::PRESSURE) {
        const _ipg_t ispgp = _espg_to_pressure.size();
        _espg_to_pressure.insert({espg, physical_quantity_value});
        _espg_to_ispg.insert({espg, ispgp});
    }
    else if (physical_quantity_type == PhysicalQuantity::VOLUME_VELOCITY) {
        const _ipg_t ispgv = _espg_to_volvel.size();
        _espg_to_volvel.insert({espg, physical_quantity_value});
        _espg_to_ispg.insert({espg, ispgv});
    }
    else {
        log::error("Possible physical quantities are volume velocity"
            " or pressure.");
    }
    _is_any_source_defined = true;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_surface_specific_acoustic_impedance(
    const _epg_t& espg,
    const _FuncRealToCmplx& impedance
) {
    _check_if_mesh_is_defined();
    if (!_existing_espg.contains(espg)) {
        log::error("Tag {} not found in mesh file.", espg);
    }
    if (_espg_to_pressure.contains(espg) ||
        _espg_to_volvel.contains(espg) ||
        _espg_to_impedance.contains(espg)
    ) {
        log::error("Tag {} already assigned.", espg);
    }
    const _ipg_t ispgi = _espg_to_impedance.size();
    _espg_to_impedance.insert({espg, impedance});
    _espg_to_ispg.insert({espg, ispgi});
}

} // namespace numav
