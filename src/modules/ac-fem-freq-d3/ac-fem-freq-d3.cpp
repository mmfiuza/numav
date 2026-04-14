// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/impl.hpp"

namespace numav {

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::set_maximum_frequency(
    const double& freq_max
) {
    pimpl->set_maximum_frequency(freq_max);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::set_frequency_range(
    const double& freq_min,
    const double& freq_max
) {
    pimpl->set_frequency_range(freq_min, freq_max);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::set_frequency_steps_count(
    const size_t& frequency_steps_count
) {
    pimpl->set_frequency_steps_count(frequency_steps_count);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::set_frequency_sampling_density(
    const FrequencySamplingDensity& frequency_sampling_density
) {
    pimpl->set_frequency_sampling_density(frequency_sampling_density);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::set_frequency_steps(
    const std::vector<double>& frequency_steps
) {
    pimpl->set_frequency_steps(frequency_steps);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::load_mesh(
    const char* const path_to_mesh
) {
    pimpl->load_mesh(path_to_mesh);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_volume_material(
    const size_t& evpg,
    const _FuncRealToCmplx& density,
    const _FuncRealToCmplx& soundspeed
) {
    pimpl->add_volume_material(evpg, density, soundspeed);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_volume_material(
    const size_t& evpg,
    const char* const density_text_file,
    const _FuncRealToCmplx& soundspeed
) {
    pimpl->add_volume_material(evpg, density_text_file, soundspeed);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_volume_material(
    const size_t& evpg,
    const _FuncRealToCmplx& density,
    const char* const soundspeed_text_file
) {
    pimpl->add_volume_material(evpg, density, soundspeed_text_file);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_volume_material(
    const size_t& evpg,
    const char* const density_text_file,
    const char* const soundspeed_text_file
) {
    pimpl->add_volume_material(evpg, density_text_file, soundspeed_text_file);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_sound_source(
    const TypeOfSource& type_of_source,
    const std::array<double,3>& point_coordinates,
    const PhysicalQuantity& physical_quantity_type,
    const _FuncRealToCmplx& physical_quantity_value
) {
    pimpl->add_sound_source(
        type_of_source,
        point_coordinates,
        physical_quantity_type,
        physical_quantity_value
    );
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_sound_source(
    const TypeOfSource& type_of_source,
    const std::array<double,3>& point_coordinates,
    const PhysicalQuantity& physical_quantity_type,
    const char* const physical_quantity_value_text_file
) {
    pimpl->add_sound_source(
        type_of_source,
        point_coordinates,
        physical_quantity_type,
        physical_quantity_value_text_file
    );
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_sound_source(
    const TypeOfSource& type_of_source,
    const size_t& espg,
    const PhysicalQuantity& physical_quantity_type,
    const _FuncRealToCmplx& physical_quantity_value
) {
    pimpl->add_sound_source(
        type_of_source,
        espg,
        physical_quantity_type,
        physical_quantity_value
    );
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_sound_source(
    const TypeOfSource& type_of_source,
    const size_t& espg,
    const PhysicalQuantity& physical_quantity_type,
    const char* const physical_quantity_value_text_file
) {
    pimpl->add_sound_source(
        type_of_source,
        espg,
        physical_quantity_type,
        physical_quantity_value_text_file
    );
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_receiver(
    const std::array<double,3>& point_coordinates
) {
    pimpl->add_receiver(point_coordinates);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_surface_specific_acoustic_impedance(
    const size_t& espg,
    const _FuncRealToCmplx& impedance
) {
    pimpl->add_surface_specific_acoustic_impedance(espg, impedance);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_surface_specific_acoustic_impedance(
    const size_t& espg,
    const char* const impedance_text_file
) {
    pimpl->add_surface_specific_acoustic_impedance(espg, impedance_text_file);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::run() {
    pimpl->run();
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::export_result(
    const char* const path_to_result
) {
    pimpl->export_result(path_to_result);
}

// explicit instantiation declarations
INSTANTIATE_SIMULATION_CLASS

} // namespace numav
