// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav.hpp"

std::complex<double> get_complex_volume_velocity_amplitude(const double& frequency_in_hz) {
    return 1/frequency_in_hz;
}

std::complex<double> get_complex_pressure_amplitude(const double& frequency_in_hz) {
    return 1;
}

std::complex<double> get_specific_surface_acoustic_impedance(const double& frequency_in_hz) {
    return 1;
}

int main() {

    using namespace numav;

    // create the simulation object with some numerical method
    auto s = Simulation<
        Phenomenon::ACOUSTIC,
        NumericalMethod::FEM,
        Domain::FREQUENCY,
        Dimension::D3
    >();
    
    // set element order
    s.set_element_order(2);

    // determine simulation frequency range
    double freq_min = 0;
    double freq_max = 100;
    s.set_freq_limits(freq_min, freq_max);

    // load the mesh
    s.load_mesh("example1.bdf");

    // add a volume material (air in this case)
    const uint64_t id_air = 1;
    const double rho_air = 1.20;
    const double c_air = 343;
    s.add_volume_material(id_air, rho_air, c_air);

    // add volume velocity sources
    std::array<double,3> source_coordinates_1 = {1.0, 1.5, 2.0};
    s.add_source(
        TypeOfSource::POINT, source_coordinates_1,
        PhysicalQuantity::VOLUME_VELOCITY, get_complex_volume_velocity_amplitude
    );
    const uint64_t id_surface_source_1 = 2;
    s.add_source(
        TypeOfSource::SURFACE, id_surface_source_1,
        PhysicalQuantity::VOLUME_VELOCITY, get_complex_volume_velocity_amplitude
    );

    // add pressure sources
    std::array<double,3> source_coordinates_2 = {2.0, 2.5, 1.0};
    s.add_source(
        TypeOfSource::POINT, source_coordinates_2,
        PhysicalQuantity::PRESSURE, get_complex_pressure_amplitude
    );
    const uint64_t id_surface_source_2 = 3;
    s.add_source(
        TypeOfSource::SURFACE, id_surface_source_2,
        PhysicalQuantity::PRESSURE, get_complex_pressure_amplitude
    );

    // add specific surface impedance
    const uint64_t id_surface_impedance = 4;
    s.add_specific_surface_acoustic_impedance(
        id_surface_impedance, get_specific_surface_acoustic_impedance
    );

    // run the simulation
    std::vector<std::vector<std::complex<double>>> complex_pressure_amplitude = s.run();
}
