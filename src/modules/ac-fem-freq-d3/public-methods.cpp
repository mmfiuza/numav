// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/impl.hpp"

#include "common/exception.hpp"
#include "common/log.hpp"
#include "common/utils.hpp"

#include <tuple>

namespace numav {

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_check_if_mesh_is_defined() {
    if (!_is_mesh_defined){
        error("Mesh not defined. Call load_mesh to do so.");
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::set_maximum_frequency(
    const Float freq_max
) {
    if (freq_max < 0_F) {
        error("Maximum frequency should be positive.");
    }
    if (freq_max == 0_F) {
        error("Maximum frequency should not be zero.");
    }
    if (_freq_type_defined_by_user != _FreqTypeDefinedByUser::UNDEFINED) {
        error("Simulation frequency is already defined.");
    }
    _freq_max = freq_max;
    _freq_type_defined_by_user = _FreqTypeDefinedByUser::MAXIMUM;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::set_frequency_range(
    const Float freq_min,
    const Float freq_max
) {
    if (freq_min < 0_F || freq_max < 0_F) {
        error("Frequency limits should be positive.");
    }
    if (freq_min == 0_F || freq_max == 0_F) {
        error("Frequency limits should not be zero.");
    }
    if (freq_min >= freq_max) {
        error("Upper frequency should be greater than the lower.");
    }
    if (_freq_type_defined_by_user != _FreqTypeDefinedByUser::UNDEFINED) {
        error("Simulation frequency is already defined.");
    }
    _freq_min = freq_min;
    _freq_max = freq_max;
    _freq_type_defined_by_user = _FreqTypeDefinedByUser::RANGE;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::set_frequency_steps_count(
    const size_t freq_steps_count
) {
    if (_freq_type_defined_by_user == _FreqTypeDefinedByUser::UNDEFINED) {
        error("Simulation frequency was not defined."
            " Call set_maximum_frequency to do so.");
    }
    if (_freq_type_defined_by_user == _FreqTypeDefinedByUser::STEPS) {
        error("Simulation frequency steps already defined.");
    }
    _freq_count = freq_steps_count;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::set_frequency_sampling_density(
    const FrequencySamplingDensity freq_sampling_density
) {
    if (_freq_type_defined_by_user == _FreqTypeDefinedByUser::UNDEFINED) {
        error("Simulation frequency was not defined."
            " Call set_maximum_frequency to do so.");
    }
    if (_freq_type_defined_by_user == _FreqTypeDefinedByUser::STEPS) {
        error("Simulation frequency steps already defined.");
    }
    _freq_sampling_density = freq_sampling_density;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::set_frequency_steps(
    const std::vector<Float>& freq_steps
) {
    if (_freq_type_defined_by_user != _FreqTypeDefinedByUser::UNDEFINED) {
        error("Simulation frequency is already defined.");
    }
    _freq_count = freq_steps.size();
    _freq_min = freq_steps.front();
    _freq_max = freq_steps.back();
    _freq_steps = fz::SafePtr<Float>(_freq_count);
    for (size_t i = 0UL; i != _freq_count; ++i) {
        _freq_steps[i] = freq_steps[i];
    }
    _freq_type_defined_by_user = _FreqTypeDefinedByUser::STEPS;

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

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::add_volume_material(
    const size_t evpg,
    const FuncFloatToCmplx& density_func,
    const FuncFloatToCmplx& soundspeed_func
) {
    _check_if_mesh_is_defined();
    if (!_existing_evpg.contains(evpg)) {
        error("Tag {} not found in mesh file.", evpg);
    }
    if (_evpg_to_volprop.contains(evpg)) {
        error("Tag {} already assigned.", evpg);
    }
    _evpg_to_volprop.insert({evpg, {density_func, soundspeed_func}});
    const size_t ivpg = _evpg_to_ivpg.size();
    _evpg_to_ivpg.insert({evpg, ivpg});
    ++_ivpg_count;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::add_volume_material(
    const size_t evpg,
    const char* const density_table,
    const FuncFloatToCmplx& soundspeed_func
) {
    FuncFloatToCmplx density_func = 
        convert_table_to_real_to_cmplx_func(density_table);

    add_volume_material(evpg, density_func, soundspeed_func);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::add_volume_material(
    const size_t evpg,
    const FuncFloatToCmplx& density_func,
    const char* const soundspeed_table
) {
    FuncFloatToCmplx soundspeed_func = 
        convert_table_to_real_to_cmplx_func(soundspeed_table);

    add_volume_material(evpg, density_func, soundspeed_func);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::add_volume_material(
    const size_t evpg,
    const char* const density_table,
    const char* const soundspeed_table
) {
    FuncFloatToCmplx density_func = 
        convert_table_to_real_to_cmplx_func(density_table);
    FuncFloatToCmplx soundspeed_func = 
        convert_table_to_real_to_cmplx_func(soundspeed_table);

    add_volume_material(evpg, density_func, soundspeed_func);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::add_sound_source(
    const TypeOfSource type_of_source,
    const std::array<Float,3UL> point_coordinates,
    const PhysicalQuantity pq_type,
    const FuncFloatToCmplx& pq_func
) {
    _check_if_mesh_is_defined();
    // TODO: check if the point is outside the mesh
    if (type_of_source != TypeOfSource::POINT) {
        error("Tried to assign coordinates to a surface sound source.");
    }
    const size_t closest_ni = _get_closest_point(point_coordinates);
    
    if (pq_type == PhysicalQuantity::VOLUME_VELOCITY) {
        _point_volvel.push_back(
            std::make_tuple(closest_ni, pq_func)
        );
        ++_pvni_count;
    }
    else if (pq_type == PhysicalQuantity::PRESSURE) {
        _point_pressure.push_back(
            std::make_tuple(closest_ni, pq_func)
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
    const TypeOfSource type_of_source,
    const std::array<Float,3UL> point_coordinates,
    const PhysicalQuantity pq_type,
    const char* const pq_table
) {
    FuncFloatToCmplx pq_func = convert_table_to_real_to_cmplx_func(pq_table);

    add_sound_source(
        type_of_source,
        point_coordinates,
        pq_type,
        pq_func
    );
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_validate_espg(const size_t espg) {
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
    const TypeOfSource type_of_source,
    const size_t espg,
    const PhysicalQuantity pq_type,
    const FuncFloatToCmplx& pq_func
) {
    _check_if_mesh_is_defined();
    _validate_espg(espg);
    if (type_of_source != TypeOfSource::SURFACE) {
        error("Tried to assign a tag to a point.");
    }
    
    if (pq_type == PhysicalQuantity::PARTICLE_VELOCITY) {
        const size_t ispgv = _espg_to_velocity.size();
        _espg_to_velocity.insert({espg, pq_func});
        _espg_to_ispg.insert({espg, ispgv});
        ++_ispgv_count;
    }
    else if (pq_type == PhysicalQuantity::PRESSURE) {
        const size_t ispgp = _espg_to_pressure.size();
        _espg_to_pressure.insert({espg, pq_func});
        _espg_to_ispg.insert({espg, ispgp});
        ++_ispgp_count;
    }
    else {
        error("Possible physical quantity types are particle velocity, "
            "volume velocity or pressure.");
    }
    _is_any_source_defined = true;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::add_sound_source(
    const TypeOfSource type_of_source,
    const size_t espg,
    const PhysicalQuantity pq_type,
    const char* const pq_table
) {
    FuncFloatToCmplx pq_func = convert_table_to_real_to_cmplx_func(pq_table);

    add_sound_source(
        type_of_source,
        espg,
        pq_type,
        pq_func
    );
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::add_receiver(
    const std::array<Float,DIM> point_coordinates
) {
    _check_if_mesh_is_defined();
    // TODO: check if the point is outside the mesh
    const Float& x = point_coordinates[0UL];
    const Float& y = point_coordinates[1UL];
    const Float& z = point_coordinates[2UL];
    _receiver_points.emplace_back(std::array<Float,DIM>{x, y, z});
    ++_ri_count;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::add_surface_material(
    const size_t espg,
    const PhysicalQuantity pq_type,
    const FuncFloatToCmplx& pq_func
) {
    _check_if_mesh_is_defined();
    _validate_espg(espg);
    if (pq_type != PhysicalQuantity::IMPEDANCE) {
        error("Possible physical quantity type is only impedance");
    }
    const size_t ispgi = _espg_to_impedance.size();
    _espg_to_impedance.insert({espg, pq_func});
    _espg_to_ispg.insert({espg, ispgi});
    ++_ispgi_count;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::add_surface_material(
    const size_t espg,
    const PhysicalQuantity pq_type,
    const char* const pq_table
) {
    FuncFloatToCmplx pq_func = convert_table_to_real_to_cmplx_func(pq_table);
    
    add_surface_material(
        espg,
        pq_type,
        pq_func
    );
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_check_if_it_can_run() {
    _check_if_mesh_is_defined();
    if (!_is_any_source_defined){
        error("No sound source was defined."
            " Call add_sound_source to do so.");
    }
    if (_freq_type_defined_by_user == _FreqTypeDefinedByUser::UNDEFINED){
        error("Simulation frequency was not defined."
            " Call set_maximum_frequency to do so.");
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
void SimulationAcFemFreqD3<O>::Impl::set_result_export_path(
    const char* const path_to_result
) {
    if (!_nmvr_file_path.empty()) {
        error("Result export path is already defined.");
    }
    _nmvr_file_path = path_to_result;
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
    _post_process();
}

// explicit instantiation declarations
INSTANTIATE_SIMULATION_CLASS

} // namespace numav
