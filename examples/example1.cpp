// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav.hpp"

int main() {

    using namespace numav;

    // create the simulation object with some numerical method
    auto s = Simulation<
        Phenomenon::ACOUSTIC,
        NumericalMethod::FEM,
        Domain::FREQUENCY,
        Dimension::D3,
        ElementOrder::O2
    >();

    // load the mesh
    s.load_mesh("./build/examples_bin/mesh1.bdff");

    // determine simulation frequency range
    double freq_min = 0;
    double freq_max = 100;
    s.set_freq_range(freq_min, freq_max);

    auto rho_air = [](auto f) { return 1.20; };
    auto c_air = [](auto f) { return 343; };
    s.add_volume_material(
        1, // id in mesh file
        rho_air,
        c_air
    );

    // add volume velocity sources
    auto q =[](auto f) { return 1/f; };
    s.add_source(
        TypeOfSource::POINT, {1.0, 1.5, 2.0}, // x,y,z
        PhysicalQuantity::VOLUME_VELOCITY, q
    );
    s.add_source(
        TypeOfSource::SURFACE, 2, // id in mesh file
        PhysicalQuantity::VOLUME_VELOCITY, q
    );

    // add pressure sources
    auto p =[](auto f) { return 1/f; };
    s.add_source(
        TypeOfSource::POINT, {2.0, 2.5, 1.0}, // x,y,z
        PhysicalQuantity::PRESSURE, p
    );
    const size_t id_surface_source_2 = 3;
    s.add_source(
        TypeOfSource::SURFACE, 3, // id in mesh file
        PhysicalQuantity::PRESSURE, p
    );

    // add specific surface acoustic impedance
    auto Z = [](auto f) { return std::complex<double>(f,2); };
    s.add_surface_specific_acoustic_impedance(
        4, // id in mesh file
        Z
    );

    // run the simulation
    auto result = s.run();
}
