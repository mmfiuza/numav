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

    template<Dimension D> constexpr size_t DIMENSION_COUNT = []{
        if constexpr (D == Dimension::D1) return 1;
        if constexpr (D == Dimension::D2) return 2;
        if constexpr (D == Dimension::D3) return 3;
        return 0; // error
    }();

    enum class ElementOrder {
        O1,
        O2
    };

    template<ElementOrder O> constexpr size_t NODES_PER_SURF_ELEMENT = []{
        if constexpr (O == ElementOrder::O1) return 3;
        if constexpr (O == ElementOrder::O2) return 6;
        return 0; // TODO: throw error here
    }();

    template<ElementOrder O> constexpr size_t NODES_PER_VOL_ELEMENT = []{
        if constexpr (O == ElementOrder::O1) return 4;
        if constexpr (O == ElementOrder::O2) return 10;
        return 0; // TODO: throw error here
    }();

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
        Phenomenon PHENOMENON, NumericalMethod NUMERICAL_METHOD, Domain DOMAIN, Dimension DIMENSION
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
        ~Result();
    private:
        double _data;
    };

    // declare the general Simulation class
    template<
        Phenomenon PHENOMENON, NumericalMethod NUMERICAL_METHOD, Domain DOMAIN, Dimension DIMENSION,
        auto... EXTRA
    >
    class Simulation {};

    // declare specific Simulation types
    template <ElementOrder ORDER>
    class Simulation<
        Phenomenon::ACOUSTIC,
        NumericalMethod::FEM,
        Domain::FREQUENCY,
        Dimension::D3,
        ORDER
    > {
    public:
        Simulation();
        ~Simulation();
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
        class Mesh {
            private:
                // TODO: remove std::vector from here
                std::vector<std::array<double,DIMENSION_COUNT<Dimension::D3>>> _node_coords;
                std::vector<std::array<size_t,NODES_PER_SURF_ELEMENT<ORDER>>> _2d_elem_vtx_idx;
                std::vector<std::array<size_t,NODES_PER_VOL_ELEMENT<ORDER>>> _3d_elem_vtx_idx;
                std::vector<size_t> _2d_elem_tag;
                std::vector<size_t> _3d_elem_tag;
                const size_t _nodes_count;
                const size_t _2d_elem_count;
                const size_t _3d_elem_count;
        };

        Mesh _mesh;
        bool _is_mesh_defined;
        bool _is_any_source_defined;
        double _freq_min;
        double _freq_max;
        std::vector<double> _freq_vector;
        std::vector<std::vector<std::complex<double>>> _complex_speed_of_sound;
        std::vector<std::vector<std::complex<double>>> _complex_density;
        std::vector<std::vector<std::complex<double>>> _specific_acoustic_impedance;
    };
} // namespace numav