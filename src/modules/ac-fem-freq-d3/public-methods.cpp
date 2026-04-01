// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/impl.hpp"

#include "common/exception.hpp"
#include "common/log.hpp"
#include "common/maths.hpp"

#include <tuple>
#include <fstream>

namespace numav {

    template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_check_if_mesh_is_defined() {
    if (!_is_mesh_defined){
        error("Mesh not defined. Call load_mesh to do so.");
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::set_frequency_range(
    const double& freq_min,
    const double& freq_max
) {
    if (freq_min<0 || freq_max<0) {
        error("Frequency limits should be positive.");
    }
    if (freq_min >= freq_max) {
        error("Upper frequency should be greater than the lower.");
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
        error("Mesh is already defined.");
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
        error("Unrecognized file format: \"{0}\".", format);
    }
    _generate_extra_nodes(); // call is based on the element order
    _is_mesh_defined = true;
}

_FuncRealToCmplx convert_table_to_real_to_cmplx_func(
    const char* const impedance_text_file
) {
    // open file
    std::ifstream file(impedance_text_file);
    std::string line;
    if (!file.is_open()) {
        error("Could not open file: {}", impedance_text_file);
    }

    // first pass: count lines
    size_t line_count = 0;
    while (std::getline(file, line)) {
        ++line_count;
    }
    file.clear();
    file.seekg(0, std::ios::beg);
    std::vector<double> real_vec;
    real_vec.reserve(line_count);
    std::vector<_cmplx_t> cmplx_vec;
    cmplx_vec.reserve(line_count);

    // second pass: read each line
    while (std::getline(file, line))
    {
        // read frequency
        size_t first_comma_pos = line.find(',');
        if (first_comma_pos == std::string::npos) {
            continue; // malformed line is skipped
        }
        std::string col1_str = line.substr(0, first_comma_pos);
        std::istringstream col1_input_string(col1_str);
        double col1;
        col1_input_string >> col1;
        real_vec.push_back(col1);
        
        // read real part of complex vector
        size_t second_comma_pos = line.find(',', first_comma_pos+1);
        if (second_comma_pos == std::string::npos) {
            continue; // malformed line is skipped
        }
        std::string col2_str =
            line.substr(first_comma_pos + 1, second_comma_pos);
        std::istringstream col2_input_string(col2_str);
        double col2;
        col2_input_string >> col2;
        
        // read imaginary part of complex vector
        std::string col3_str = line.substr(second_comma_pos + 1);
        std::istringstream col3_input_string(col3_str);
        double col3;
        col3_input_string >> col3;
        
        // write complex vector
        cmplx_vec.push_back(_cmplx_t(col2, col3));
    }
    
    // create the _FuncRealToCmplx funciton
    auto func_real_to_cmplx = [real_vec, cmplx_vec](double real) {
        return interpolate(real_vec, cmplx_vec, real);
    };

    return func_real_to_cmplx;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::add_volume_material(
    const size_t& evpg,
    const _FuncRealToCmplx& density,
    const _FuncRealToCmplx& soundspeed
) {
    _check_if_mesh_is_defined();
    if (!_existing_evpg.contains(evpg)) {
        error("Tag {} not found in mesh file.", evpg);
    }
    if (_evpg_to_volprop.contains(evpg)) {
        error("Tag {} already assigned.", evpg);
    }
    _evpg_to_volprop.insert({evpg, {density, soundspeed}});
    const size_t ivpg = _evpg_to_ivpg.size();
    _evpg_to_ivpg.insert({evpg, ivpg});
    ++_ivpg_count;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::add_volume_material(
    const size_t& evpg,
    const char* const density_text_file,
    const _FuncRealToCmplx& soundspeed
) {
    _FuncRealToCmplx density = 
        convert_table_to_real_to_cmplx_func(density_text_file);

    add_volume_material(evpg, density, soundspeed);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::add_volume_material(
    const size_t& evpg,
    const _FuncRealToCmplx& density,
    const char* const soundspeed_text_file
) {
    _FuncRealToCmplx soundspeed = 
        convert_table_to_real_to_cmplx_func(soundspeed_text_file);

    add_volume_material(evpg, density, soundspeed);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::add_volume_material(
    const size_t& evpg,
    const char* const density_text_file,
    const char* const soundspeed_text_file
) {
    _FuncRealToCmplx density = 
        convert_table_to_real_to_cmplx_func(density_text_file);
    _FuncRealToCmplx soundspeed = 
        convert_table_to_real_to_cmplx_func(soundspeed_text_file);

    add_volume_material(evpg, density, soundspeed);
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
        error("Tried to assign coordinates to a surface sound source.");
    }
    const size_t closest_ni = _get_closest_point(point_coordinates);
    
    if (physical_quantity_type == PhysicalQuantity::VOLUME_VELOCITY) {
        _point_volvel.push_back(
            std::make_tuple(closest_ni, physical_quantity_value)
        );
        ++_pvni_count;
    }
    else if (physical_quantity_type == PhysicalQuantity::PRESSURE) {
        _point_pressure.push_back(
            std::make_tuple(closest_ni, physical_quantity_value)
        );
        ++_ppni_count;
    }
    else {
        error("Possible physical quantities are volume velocity"
            " or pressure.");
    }
    _is_any_source_defined = true;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::add_sound_source(
    const TypeOfSource& type_of_source,
    const std::array<double,3>& point_coordinates,
    const PhysicalQuantity& physical_quantity_type,
    const char* const physical_quantity_value_text_file
) {
    _FuncRealToCmplx physical_quantity_value = 
        convert_table_to_real_to_cmplx_func(physical_quantity_value_text_file);

    add_sound_source(
        type_of_source,
        point_coordinates,
        physical_quantity_type,
        physical_quantity_value
    );
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_validate_espg(const size_t& espg) {
    if (!_existing_espg.contains(espg)) {
        error("Tag {} not found in mesh file.", espg);
    }
    if (_espg_to_pressure.contains(espg) ||
        _espg_to_velocity.contains(espg) ||
        _espg_to_impedance.contains(espg)
    ) {
        error("Tag {} already assigned.", espg);
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::add_sound_source(
    const TypeOfSource& type_of_source,
    const size_t& espg,
    const PhysicalQuantity& physical_quantity_type,
    const _FuncRealToCmplx& physical_quantity_value
) {
    _check_if_mesh_is_defined();
    _validate_espg(espg);
    if (type_of_source != TypeOfSource::SURFACE) {
        error("Tried to assign a tag to a point.");
    }
    
    if (physical_quantity_type == PhysicalQuantity::PARTICLE_VELOCITY) {
        const size_t ispgv = _espg_to_velocity.size();
        _espg_to_velocity.insert({espg, physical_quantity_value});
        _espg_to_ispg.insert({espg, ispgv});
        ++_ispgv_count;
    }
    else if (physical_quantity_type == PhysicalQuantity::PRESSURE) {
        const size_t ispgp = _espg_to_pressure.size();
        _espg_to_pressure.insert({espg, physical_quantity_value});
        _espg_to_ispg.insert({espg, ispgp});
        ++_ispgp_count;
    }
    else {
        error("Possible physical quantities are particle velocity"
            " or pressure.");
    }
    _is_any_source_defined = true;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::add_sound_source(
    const TypeOfSource& type_of_source,
    const size_t& espg,
    const PhysicalQuantity& physical_quantity_type,
    const char* const physical_quantity_value_text_file
) {
    _FuncRealToCmplx physical_quantity_value = 
        convert_table_to_real_to_cmplx_func(physical_quantity_value_text_file);

    add_sound_source(
        type_of_source,
        espg,
        physical_quantity_type,
        physical_quantity_value
    );
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::add_surface_specific_acoustic_impedance(
    const size_t& espg,
    const _FuncRealToCmplx& impedance
) {
    _check_if_mesh_is_defined();
    _validate_espg(espg);
    const size_t ispgi = _espg_to_impedance.size();
    _espg_to_impedance.insert({espg, impedance});
    _espg_to_ispg.insert({espg, ispgi});
    ++_ispgi_count;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::add_surface_specific_acoustic_impedance(
    const size_t& espg,
    const char* const impedance_text_file
) {
    _FuncRealToCmplx impedance = 
        convert_table_to_real_to_cmplx_func(impedance_text_file);
    
    add_surface_specific_acoustic_impedance(espg, impedance);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_check_if_it_can_run() {
    _check_if_mesh_is_defined();
    if (!_is_any_source_defined){
        error("No sound source was defined."
            " Call add_sound_source to do so.");
    }
    if (!_is_freq_range_defined){
        error("Simulation frequency range was not defined."
            " Call set_frequency_range to do so.");
    }
    for (auto& evpg : _existing_evpg) {
        if (!_evpg_to_volprop.contains(evpg)) {
            error("Volume tag {} was not assigned."
            " Call add_volume_material to do so.", evpg);
        }
    }
    if (_did_run) {
        error("This Simulation has already been run.");
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::run()
{
    _check_if_it_can_run();
    log::print_opening();
    log::print_opening_ac_fem_freq_d3();
    _define_freq_vector();
    _organize_physical_group_data();
    _assemble_freq_independent_parts();
    _solve_systems();
}

// TODO: move export_result to here

// explicit instantiation declarations
INSTANTIATE_SIMULATION_CLASS

} // namespace numav
