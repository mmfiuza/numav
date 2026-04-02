// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav/numav.hpp"

int main()
{
    using namespace numav;

    // create the simulation object with some numerical method
    auto s = Simulation<
        Phenomenon::ACOUSTIC,
        NumericalMethod::FEM,
        Domain::FREQUENCY,
        Dimension::D3,
        ElementOrder::O1
    >();

    // load the mesh
    s.load_mesh("build/tests_bin/test1.bdf");

    // determine simulation frequency range
    double freq_min = 20;
    double freq_max = 100;
    s.set_frequency_range(freq_min, freq_max);

    s.add_volume_material(
        1,
        [](auto f) { return 1.2; },
        [](auto f) { return std::complex<double>(343,0); }
    );
    s.add_volume_material(
        2,
        [](auto f) { return 1.2; },
        [](auto f) { return std::complex<double>(343,5); }
    );

    // add volume velocity sources
    auto q = [](auto f) { return 1/f; };
    s.add_sound_source(
        TypeOfSource::POINT, {3.0, 2.0, 1.0},
        PhysicalQuantity::VOLUME_VELOCITY, q
    );
    auto u = [](auto f) { return 1/f; };
    s.add_sound_source(
        TypeOfSource::SURFACE, 2,
        PhysicalQuantity::PARTICLE_VELOCITY, u
    );

    // add pressure sources
    auto p0 = [](auto f) { return 1; };
    s.add_sound_source(
        TypeOfSource::POINT, {4.0, 1.0, 0.5},
        PhysicalQuantity::PRESSURE, p0
    );
    auto p1 = [](auto f) { return 2; };
    s.add_sound_source(
        TypeOfSource::SURFACE, 3,
        PhysicalQuantity::PRESSURE, p1
    );
 
    // add specific surface acoustic impedance
    auto Z = [](auto f) { return std::complex<double>(4000, 4000); };
    s.add_surface_specific_acoustic_impedance(5, Z);

    // run the simulation
    s.run();

    // export the result to binary
    s.export_result("result.nmvr");
}
