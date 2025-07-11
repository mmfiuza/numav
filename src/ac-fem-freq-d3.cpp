// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav.hpp"

using ResultAcFemFreqD3 = typename numav::Result<
    numav::Phenomenon::ACOUSTIC,
    numav::NumericalMethod::FEM,
    numav::Domain::FREQUENCY,
    numav::Dimension::D3
>;

ResultAcFemFreqD3::Result() = default;
ResultAcFemFreqD3::~Result() = default;

template<numav::ElementOrder O>
using SimulationAcFemFreqD3 = typename numav::Simulation<
    numav::Phenomenon::ACOUSTIC,
    numav::NumericalMethod::FEM,
    numav::Domain::FREQUENCY,
    numav::Dimension::D3,
    O
>;

template<numav::ElementOrder O>
SimulationAcFemFreqD3<O>::Simulation::Mesh::Mesh() = default;

template<numav::ElementOrder O>
SimulationAcFemFreqD3<O>::Simulation::Mesh::~Mesh() = default;

template<numav::ElementOrder O>
SimulationAcFemFreqD3<O>::Simulation() {
}

template<numav::ElementOrder O>
SimulationAcFemFreqD3<O>::~Simulation() {
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::set_freq_limits(
    const double& freq_min, const double& freq_max
) {
    _freq_min = freq_min;
    _freq_max = freq_max;
    // TODO: genereate _freq_vector
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::load_mesh(const char* const path_to_mesh) {
    
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::add_volume_material(
    const uint64_t& id, const double& rho, const double& c
) {
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::add_source(
    const TypeOfSource& type_of_source,
    const std::array<double,3>& point_coordinates,
    const PhysicalQuantity& physical_quantity,
    const std::function<std::complex<double>(double)>& physical_quantity_value
) {
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::add_source(
    const TypeOfSource& type_of_source,
    const uint64_t& surface_id,
    const PhysicalQuantity& physical_quantity,
    const std::function<std::complex<double>(double)>& physical_quantity_value
) {
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::add_surface_specific_acoustic_impedance(
    const uint64_t& surface_id,
    const std::function<std::complex<double>(double)>& impedance
) {
}

template <numav::ElementOrder O>
ResultAcFemFreqD3 SimulationAcFemFreqD3<O>::run() {
    return ResultAcFemFreqD3();
}
