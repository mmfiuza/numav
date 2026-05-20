// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/impl.hpp"

namespace numav {

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::set_maximum_frequency(
    const Float freq_max
) {
    _pimpl->set_maximum_frequency(freq_max);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::set_frequency_range(
    const Float freq_min,
    const Float freq_max
) {
    _pimpl->set_frequency_range(freq_min, freq_max);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::set_frequency_steps_count(
    const size_t freq_steps_count
) {
    _pimpl->set_frequency_steps_count(freq_steps_count);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::set_frequency_sampling_density(
    const FrequencySamplingDensity freq_sampling_density
) {
    _pimpl->set_frequency_sampling_density(freq_sampling_density);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::set_frequency_steps(
    const std::vector<Float>& freq_steps
) {
    _pimpl->set_frequency_steps(freq_steps);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::load_mesh(
    const char* const path_to_mesh
) {
    _pimpl->load_mesh(path_to_mesh);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_volume_material(
    const size_t evpg,
    const FuncFloatToCmplx& density_func,
    const FuncFloatToCmplx& soundspeed_func
) {
    _pimpl->add_volume_material(evpg, density_func, soundspeed_func);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_volume_material(
    const size_t evpg,
    const char* const density_table,
    const FuncFloatToCmplx& soundspeed_func
) {
    _pimpl->add_volume_material(evpg, density_table, soundspeed_func);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_volume_material(
    const size_t evpg,
    const FuncFloatToCmplx& density_func,
    const char* const soundspeed_table
) {
    _pimpl->add_volume_material(evpg, density_func, soundspeed_table);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_volume_material(
    const size_t evpg,
    const char* const density_table,
    const char* const soundspeed_table
) {
    _pimpl->add_volume_material(evpg, density_table, soundspeed_table);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_sound_source(
    const TypeOfSource source_type,
    const std::array<Float,3UL> point_coords,
    const PhysicalQuantity pq_type,
    const FuncFloatToCmplx& pq_func
) {
    _pimpl->add_sound_source(source_type, point_coords, pq_type, pq_func);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_sound_source(
    const TypeOfSource source_type,
    const std::array<Float,3UL> point_coords,
    const PhysicalQuantity pq_type,
    const char* const pq_table
) {
    _pimpl->add_sound_source(source_type, point_coords, pq_type, pq_table);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_sound_source(
    const TypeOfSource source_type,
    const size_t espg,
    const PhysicalQuantity pq_type,
    const FuncFloatToCmplx& pq_func
) {
    _pimpl->add_sound_source(source_type, espg, pq_type, pq_func);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_sound_source(
    const TypeOfSource source_type,
    const size_t espg,
    const PhysicalQuantity pq_type,
    const char* const pq_table
) {
    _pimpl->add_sound_source(source_type, espg, pq_type, pq_table);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_receiver(
    const std::array<Float,3UL> point_coords
) {
    _pimpl->add_receiver(point_coords);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_surface_material(
    const size_t espg,
    const PhysicalQuantity pq_type,
    const FuncFloatToCmplx& pq_func
) {
    _pimpl->add_surface_material(espg, pq_type, pq_func);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_surface_material(
    const size_t espg,
    const PhysicalQuantity pq_type,
    const char* const pq_table
) {
    _pimpl->add_surface_material(espg, pq_type, pq_table);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::set_result_export_path(
    const char* const path_to_result
) {
    _pimpl->set_result_export_path(path_to_result);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::run() {
    _pimpl->run();
}

// explicit instantiation declarations
INSTANTIATE_SIMULATION_CLASS

} // namespace numav
