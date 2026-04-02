// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/exception.hpp"
#include "common/nmvr-format.hpp"
#include "common/utils.hpp"
#include "modules/ac-fem-freq-d3/impl.hpp"

namespace numav {

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_write_nmvr(
    const char* const path_to_result
) {
    std::ofstream file(path_to_result, std::ios::binary);
    if (!file) {
        error("Failed to open file: {}", path_to_result);
    }

    // signature
    file.write(nmvr::SIGNATURE.data(), nmvr::SIGNATURE_SIZE);

    // chunk: simulation type
    constexpr uint64_t SIMULATION_TYPE_CHUNK_SIZE = // bytes
        nmvr::SIM_TYPE_PHENOMENON_SIZE +
        nmvr::SIM_TYPE_NUMERICAL_METHOD_SIZE +
        nmvr::SIM_TYPE_DOMAIN_SIZE +
        nmvr::SIM_TYPE_DIMENSION_SIZE +
        nmvr::SIM_TYPE_ELEMENT_ORDER_SIZE;
    constexpr std::array<char, SIMULATION_TYPE_CHUNK_SIZE> 
        SIMULATION_TYPE_DATA =
            concat_constexpr_arrays(
                nmvr::SIM_TYPE_PHENOMENON,
                nmvr::SIM_TYPE_NUMERICAL_METHOD,
                nmvr::SIM_TYPE_DOMAIN<Domain::FREQUENCY>,
                nmvr::SIM_TYPE_DIMENSION<Dimension::D3>,
                nmvr::SIM_TYPE_FEM_ORDER<O>
            );
    nmvr::write_data_chunk(
        file,
        nmvr::SIMULATION_TYPE_CHUNK_ID,
        SIMULATION_TYPE_CHUNK_SIZE,
        SIMULATION_TYPE_DATA.data()
    );

    // chunk: simulated frequency steps
    nmvr::write_data_chunk(
        file,
        nmvr::SIMULATED_FREQUENCY_STEPS_CHUNK_ID,
        _freq_count * sizeof(uint64_t),
        _freq_steps.data()
    );

    // chunk: node index to coordinates
    nmvr::write_data_chunk(
        file,
        nmvr::NODE_INDEX_TO_COORDINATES_CHUNK_ID,
        _ni_count * 3 * sizeof(double),
        _ni_to_coords.data()
    );

    // chunk: surface element index to node index
    nmvr::write_data_chunk(
        file,
        nmvr::SURFACE_ELEMENT_INDEX_TO_NODE_INDEX_CHUNK_ID,
        _sei_count * NODES_IN_SFC_ELEM<O> * sizeof(uint64_t),
        _sei_to_ni.data()
    );

    // chunk: volume element index to node index
    nmvr::write_data_chunk(
        file,
        nmvr::VOLUME_ELEMENT_INDEX_TO_NODE_INDEX_CHUNK_ID,
        _vei_count * NODES_IN_VOL_ELEM<O> * sizeof(uint64_t),
        _vei_to_ni.data()
    );

    //chunk: complex amplitude of sound pressure solution
    nmvr::write_data_chunk(
        file,
        nmvr::COMPLEX_AMPLITUDE_OF_SOUND_PRESSURE_SOLUTION_CHUNK_ID,
        _ni_count * _freq_count * sizeof(_cmplx_t),
        _cmplx_pressure_amp.data()
    );

    // end of file
    file.write(nmvr::END_OF_FILE.data(), nmvr::END_OF_FILE_SIZE);

    file.close();
}

// explicit instantiation declarations
INSTANTIATE_SIMULATION_CLASS

} // namespace numav
