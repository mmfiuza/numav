// Copyright (c) 2025 Matheus Machado Fiuza <matheus.fiuza@eac.ufsm.br>

#pragma once

#include <array>
#include <vector>
#include <complex>
#include <functional>

namespace numav {

    enum class Phenomenon {
        acoustic,
    };

    enum class NumericalMethod {
        fem,
    };

    enum class Domain {
        frequency,
        time
    };

    enum class Dimension {
        d1,
        d2,
        d3,
    };

    template <Phenomenon phe, NumericalMethod met, Domain dom, Dimension dim>
    class Simulation {};

    enum TypeOfSource {
        point,
        surface,
    };

    enum PhysicalQuantity {
        pressure,
        volume_velocity,
    };

    template <> class Simulation<
        Phenomenon::acoustic,
        NumericalMethod::fem,
        Domain::frequency,
        Dimension::d3
    > {
    public:
        Simulation();
        ~Simulation();
        void set_element_order(const uint64_t&);
        void set_freq_limits(const double&, const double&);
        void load_mesh(const std::string&);
        void add_volume_material(const uint64_t&, const double&, const double&);
        void add_source(
            const TypeOfSource&, const std::array<double,3>&,
            const PhysicalQuantity&, const std::function<std::complex<double>(double)>&
        );
        void add_source(
            const TypeOfSource&, const uint64_t&,
            const PhysicalQuantity&, const std::function<std::complex<double>(double)>&
        );
        void add_specific_surface_acoustic_impedance(
            const uint64_t&, const std::function<std::complex<double>(double)>&
        );
        std::vector<std::vector<std::complex<double>>> run();
    private:
        // Mesh _mesh;
        uint64_t _freq_min;
        uint64_t _freq_max;
        std::vector<double> _freq_vector;
        std::vector<std::vector<std::complex<double>>> _complex_speed_of_sound;
        std::vector<std::vector<std::complex<double>>> _complex_density;
        std::vector<std::vector<std::complex<double>>> _specific_acoustic_impedance;
    };
}