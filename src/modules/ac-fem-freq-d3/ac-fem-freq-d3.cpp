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

#define NUMAV_ADD_SOUND_SOURCE_PIMPL(TYPE1, TYPE2) \
    template <ElementOrder O> \
    void SimulationAcFemFreqD3<O>::add_sound_source( \
        const TypeOfSource source_type, \
        TYPE1 position, \
        const PhysicalQuantity pq_type, \
        TYPE2 pq_value \
    ) { \
        _pimpl->add_sound_source(source_type, position, pq_type, pq_value); \
    }
NUMAV_ADD_SOUND_SOURCE_PIMPL(NUMAV_ARRAY3, const FuncFloatToCmplx&)
NUMAV_ADD_SOUND_SOURCE_PIMPL(NUMAV_ARRAY3, const Cmplx            )
NUMAV_ADD_SOUND_SOURCE_PIMPL(NUMAV_ARRAY3, const char* const      )
NUMAV_ADD_SOUND_SOURCE_PIMPL(const size_t, const FuncFloatToCmplx&)
NUMAV_ADD_SOUND_SOURCE_PIMPL(const size_t, const Cmplx            )
NUMAV_ADD_SOUND_SOURCE_PIMPL(const size_t, const char* const      )

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_receiver(
    const std::array<Float,3UL> point_coords
) {
    _pimpl->add_receiver(point_coords);
}

#define NUMAV_ADD_SURFACE_MATERIAL_PIMPL(TYPE) \
    template <ElementOrder O> \
    void SimulationAcFemFreqD3<O>::add_surface_material( \
        const size_t espg, \
        const PhysicalQuantity pq_type, \
        TYPE pq_value \
    ) { \
        _pimpl->add_surface_material(espg, pq_type, pq_value); \
    }
NUMAV_ADD_SURFACE_MATERIAL_PIMPL(const FuncFloatToCmplx&)
NUMAV_ADD_SURFACE_MATERIAL_PIMPL(const Cmplx            )
NUMAV_ADD_SURFACE_MATERIAL_PIMPL(const char* const      )

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
