// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav.hpp"

using ResultAcFemFreqD3 = typename numav::Result<
    numav::Phenomenon::ACOUSTIC,
    numav::NumericalMethod::FEM,
    numav::Domain::FREQUENCY,
    numav::Dimension::D3
>;

ResultAcFemFreqD3::Result() {

}

ResultAcFemFreqD3::~Result() {

}

using SimulationAcFemFreqD3 = typename numav::Simulation<
    numav::Phenomenon::ACOUSTIC,
    numav::NumericalMethod::FEM,
    numav::Domain::FREQUENCY,
    numav::Dimension::D3
>;

SimulationAcFemFreqD3::Simulation() {

}

SimulationAcFemFreqD3::~Simulation() {

}

void SimulationAcFemFreqD3::set_element_order(const uint64_t& order) {

}

void SimulationAcFemFreqD3::set_freq_limits(const double& freq_min, const double& freq_max) {

}

void SimulationAcFemFreqD3::load_mesh(const std::string& path_to_mesh) {
    
}

void SimulationAcFemFreqD3::add_volume_material(
    const uint64_t& id, const double& rho, const double& c
) {

}

void SimulationAcFemFreqD3::add_source(
    const TypeOfSource& type_of_source,
    const std::array<double,3>& point_coordinates,
    const PhysicalQuantity& physical_quantity,
    const std::function<std::complex<double>(double)>& physical_quantity_value
) {

}

void SimulationAcFemFreqD3::add_source(
    const TypeOfSource& type_of_source,
    const uint64_t& surface_id,
    const PhysicalQuantity& physical_quantity,
    const std::function<std::complex<double>(double)>& physical_quantity_value
) {

}

void SimulationAcFemFreqD3::add_surface_specific_acoustic_impedance(
    const uint64_t& surface_id, const std::function<std::complex<double>(double)>& impedance
) {

}

ResultAcFemFreqD3 SimulationAcFemFreqD3::run() {
    return ResultAcFemFreqD3();
}
