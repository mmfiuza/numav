// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include <array>
#include <vector>
#include <complex>
#include <functional>

namespace numav {

    enum class Phenomenon {
        ACOUSTIC,
    };

    enum class NumericalMethod {
        FEM,
    };

    enum class Domain {
        FREQUENCY,
        TIME
    };

    enum class Dimension {
        D1,
        D2,
        D3,
    };

    enum class TypeOfSource {
        POINT,
        SURFACE,
    };

    enum class PhysicalQuantity {
        PRESSURE,
        VOLUME_VELOCITY,
    };

    // declare the general Result class
    template<
        Phenomenon PHENOMENON,
        NumericalMethod NUMERICAL_METHOD,
        Domain DOMAIN,
        Dimension DIMENSION
    >
    class Result {};

    template <> class Result<
        Phenomenon::ACOUSTIC,
        NumericalMethod::FEM,
        Domain::FREQUENCY,
        Dimension::D3
    > {
    public:
        Result();
    private:
        double _data;
    };

    // declare the general Simulation class
    template<
        Phenomenon PHENOMENON,
        NumericalMethod NUMERICAL_METHOD,
        Domain DOMAIN,
        Dimension DIMENSION
    >
    class Simulation {};

    // declare specific Simulation types
    template <> class Simulation<
        Phenomenon::ACOUSTIC,
        NumericalMethod::FEM,
        Domain::FREQUENCY,
        Dimension::D3
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
        void add_surface_specific_acoustic_impedance(
            const uint64_t&, const std::function<std::complex<double>(double)>&
        );
        Result<Phenomenon::ACOUSTIC, NumericalMethod::FEM, Domain::FREQUENCY, Dimension::D3> run();
    private:
        // Mesh _mesh;
        uint64_t _freq_min;
        uint64_t _freq_max;
        std::vector<double> _freq_vector;
        std::vector<std::vector<std::complex<double>>> _complex_speed_of_sound;
        std::vector<std::vector<std::complex<double>>> _complex_density;
        std::vector<std::vector<std::complex<double>>> _specific_acoustic_impedance;
    };
} // namespace numav