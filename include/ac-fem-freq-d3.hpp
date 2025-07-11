// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

namespace numav {

template<> class Result<
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

template<ElementOrder O>
class Simulation<
    Phenomenon::ACOUSTIC,
    NumericalMethod::FEM,
    Domain::FREQUENCY,
    Dimension::D3,
    O
> {
public:
    Simulation();
    ~Simulation();
    void set_freq_limits(
        const double&,
        const double&
    );
    void load_mesh(
        const char* const
    );
    void add_volume_material(
        const uint64_t&,
        const double&,
        const double&
    );
    void add_source(
        const TypeOfSource&,
        const std::array<double,3>&,
        const PhysicalQuantity&,
        const std::function<std::complex<double>(double)>&
    );
    void add_source(
        const TypeOfSource&,
        const uint64_t&,
        const PhysicalQuantity&,
        const std::function<std::complex<double>(double)>&
    );
    void add_surface_specific_acoustic_impedance(
        const uint64_t&, const std::function<std::complex<double>(double)>&
    );
    Result<
        Phenomenon::ACOUSTIC,
        NumericalMethod::FEM,
        Domain::FREQUENCY,
        Dimension::D3
    > run();

private:
    class Mesh {
        public:
        Mesh();
        ~Mesh();

        private:
        SafePtr<std::array<double,DIMENSION_COUNT<Dimension::D3>>> _node_coords;
        SafePtr<std::array<size_t,NODES_PER_SURF_ELEMENT<O>>> _2d_elem_vtx_idx;
        SafePtr<std::array<size_t,NODES_PER_VOL_ELEMENT<O>>> _3d_elem_vtx_idx;
        SafePtr<size_t> _2d_elem_tag;
        SafePtr<size_t> _3d_elem_tag;
        size_t _nodes_count;
        size_t _2d_elem_count;
        size_t _3d_elem_count;
    };

    Mesh _mesh;
    bool _is_mesh_defined;
    bool _is_any_source_defined;
    double _freq_min;
    double _freq_max;

    std::vector<double> _freq_vector;

    std::vector<std::vector<std::complex<double>>>
    _complex_speed_of_sound;

    std::vector<std::vector<std::complex<double>>>
    _complex_density;

    std::vector<std::vector<std::complex<double>>>
    _specific_acoustic_impedance;
};

// explicit instantiation declarations
template class Simulation<
    Phenomenon::ACOUSTIC,
    NumericalMethod::FEM,
    Domain::FREQUENCY,
    Dimension::D3,
    ElementOrder::O1
>;
template class Simulation<
    Phenomenon::ACOUSTIC,
    NumericalMethod::FEM,
    Domain::FREQUENCY,
    Dimension::D3,
    ElementOrder::O2
>;

} // namespace numav
