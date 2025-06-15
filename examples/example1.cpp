// Copyright (c) 2025 Matheus Machado Fiuza <matheus.fiuza@eac.ufsm.br>

#include "numav.hpp"

std::complex<double> get_complex_volume_velocity_amplitude(const double& frequnecy_in_hz) {
    return 1/frequnecy_in_hz;
}

std::complex<double> get_complex_pressure_amplitude(const double& frequnecy_in_hz) {
    return 1;
}

std::complex<double> get_specific_surface_acoustic_impedance(const double& frequnecy_in_hz) {
    return 1;
}

int main() {

    using namespace numav;

    // create the simulation object with some numerical method
    auto sml = Simulation<acoustic,fem,freq_domain,d3>();
    
    // set element order
    sml.set_element_order(2);

    // determine simulation frequency range
    double freq_min = 0;
    double freq_max = 100;
    sml.set_freq_limits(freq_min, freq_max);

    // load the mesh
    sml.load_mesh("example1.bdf");

    // add a volume material (air in this case)
    uint64_t id_air = 1;
    double rho_air = 1.20;
    double c_air = 343;
    sml.add_volume_material(id_air, rho_air, c_air);

    // add volume velocity sources
    std::array<double,3> source_coordinates_1 = {1.0, 1.5, 2.0};
    sml.add_source(
        point, source_coordinates_1, volume_velocity, get_complex_volume_velocity_amplitude
    );
    uint64_t id_surface_source_1 = 2;
    sml.add_source(
        surface, id_surface_source_1, volume_velocity, get_complex_volume_velocity_amplitude
    );

    // add pressure sources
    std::array<double,3> source_coordinates_2 = {2.0, 2.5, 1.0};
    sml.add_source(
        point, source_coordinates_2, pressure, get_complex_pressure_amplitude
    );
    uint64_t id_surface_source_2 = 3;
    sml.add_source(
        surface, id_surface_source_2, pressure, get_complex_pressure_amplitude
    );

    // add specific surface impedance
    uint64_t id_surface_impedance = 4;
    sml.add_specific_surface_acoustic_impedance(
        id_surface_impedance, get_specific_surface_acoustic_impedance
    );

    // run the simulation
    std::vector<std::vector<std::complex<double>>> complex_pressure_amplitude = sml.run();
}