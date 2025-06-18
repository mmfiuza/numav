// Copyright (c) 2025 Matheus Machado Fiuza <matheus.fiuza@eac.ufsm.br>

#include "numav.hpp"

using AcFemFreqD3 = typename numav::Simulation<
    numav::Phenomenon::ACOUSTIC,
    numav::NumericalMethod::FEM,
    numav::Domain::FREQUENCY,
    numav::Dimension::D3
>;

AcFemFreqD3::Simulation() {

}

void AcFemFreqD3::set_element_order(const uint64_t& order) {

}

void AcFemFreqD3::set_freq_limits(
    const double& freq_min, const double& freq_max) {

}

void AcFemFreqD3::load_mesh(const std::string& path_to_mesh) {
    
}

void AcFemFreqD3::add_volume_material(const uint64_t& id, const double& rho, const double& c) {

}

void AcFemFreqD3::add_source(
    const TypeOfSource& type_of_source,
    const std::array<double,3>& point_coordinates,
    const PhysicalQuantity& physical_quantity,
    const std::function<double(double)>& physical_quantity_value
) {

}
