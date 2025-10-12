// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav/numav.hpp"
#include "modules/ac-fem-freq-d3/simulation/implementation.hpp"

#include "common/log.hpp"

#include <tuple>
#include <fstream>

namespace numav {

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::set_frequency_range(
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
void SimulationAcFemFreqD3<O>::Impl::load_mesh(
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
void SimulationAcFemFreqD3<O>::Impl::add_volume_material(
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
void SimulationAcFemFreqD3<O>::Impl::add_sound_source(
    const TypeOfSource& type_of_source,
    const std::array<double,3>& point_coordinates,
    const PhysicalQuantity& physical_quantity_type,
    const _FuncRealToCmplx& physical_quantity_value
) {
    _check_if_mesh_is_defined();
    if (type_of_source != TypeOfSource::POINT) {
        log::error("Tried to assign coordinates to a surface sound source.");
    }
    const _idx_t closest_ni = _get_closest_point(point_coordinates);
    
    if (physical_quantity_type == PhysicalQuantity::VOLUME_VELOCITY) {
        _point_volvel.push_back(
            std::make_tuple(closest_ni, physical_quantity_value)
        );
    }
    else if (physical_quantity_type == PhysicalQuantity::PRESSURE) {
        _point_pressure.push_back(
            std::make_tuple(closest_ni, physical_quantity_value)
        );
    }
    else {
        log::error("Possible physical quantities are volume velocity"
            " or pressure.");
    }
    _is_any_source_defined = true;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::add_sound_source(
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
        _espg_to_velocity.contains(espg) ||
        _espg_to_impedance.contains(espg)
    ) {
        log::error("Tag {} already assigned.", espg);
    }
    if (physical_quantity_type == PhysicalQuantity::PARTICLE_VELOCITY) {
        const _ipg_t ispgv = _espg_to_velocity.size();
        _espg_to_velocity.insert({espg, physical_quantity_value});
        _espg_to_ispg.insert({espg, ispgv});
    }
    else if (physical_quantity_type == PhysicalQuantity::PRESSURE) {
        const _ipg_t ispgp = _espg_to_pressure.size();
        _espg_to_pressure.insert({espg, physical_quantity_value});
        _espg_to_ispg.insert({espg, ispgp});
    }
    else {
        log::error("Possible physical quantities are particle velocity"
            " or pressure.");
    }
    _is_any_source_defined = true;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::add_surface_specific_acoustic_impedance(
    const _epg_t& espg,
    const _FuncRealToCmplx& impedance
) {
    _check_if_mesh_is_defined();
    if (!_existing_espg.contains(espg)) {
        log::error("Tag {} not found in mesh file.", espg);
    }
    if (_espg_to_pressure.contains(espg) ||
        _espg_to_velocity.contains(espg) ||
        _espg_to_impedance.contains(espg)
    ) {
        log::error("Tag {} already assigned.", espg);
    }
    const _ipg_t ispgi = _espg_to_impedance.size();
    _espg_to_impedance.insert({espg, impedance});
    _espg_to_ispg.insert({espg, ispgi});
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::set_frequency_range(
    const double& freq_min,
    const double& freq_max
) {
    pimpl->set_frequency_range(freq_min, freq_max);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::load_mesh(
    const char* const path_to_mesh
) {
    pimpl->load_mesh(path_to_mesh);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_volume_material(
    const _epg_t& evpg,
    const _FuncRealToCmplx& density,
    const _FuncRealToCmplx& soundspeed
) {
    pimpl->add_volume_material(evpg, density, soundspeed);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_sound_source(
    const TypeOfSource& type_of_source,
    const std::array<double,3>& point_coordinates,
    const PhysicalQuantity& physical_quantity_type,
    const _FuncRealToCmplx& physical_quantity_value
) {
    pimpl->add_sound_source(
        type_of_source, point_coordinates,
        physical_quantity_type, physical_quantity_value
    );
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_sound_source(
    const TypeOfSource& type_of_source,
    const _epg_t& espg,
    const PhysicalQuantity& physical_quantity_type,
    const _FuncRealToCmplx& physical_quantity_value
) {
    pimpl->add_sound_source(
        type_of_source, espg,
        physical_quantity_type, physical_quantity_value
    );
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_surface_specific_acoustic_impedance(
    const _epg_t& espg,
    const _FuncRealToCmplx& impedance
) {
    pimpl->add_surface_specific_acoustic_impedance(espg, impedance);
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
