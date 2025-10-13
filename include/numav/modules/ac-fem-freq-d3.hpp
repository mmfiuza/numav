// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include <memory>

#include "Eigen/Eigen"

namespace numav {

template<ElementOrder O> constexpr size_t NODES_IN_SFC_ELEM = [] {
    if constexpr (O == ElementOrder::O1) { return 3; }
    if constexpr (O == ElementOrder::O2) { return 6; }
    return 0;
}();

template<ElementOrder O> constexpr size_t EXTRA_NODES_IN_SFC_ELEM = [] {
    return NODES_IN_SFC_ELEM<O> - NODES_IN_SFC_ELEM<ElementOrder::O1>;
}();

template<ElementOrder O> constexpr size_t NODES_IN_VOL_ELEM = [] {
    if constexpr (O == ElementOrder::O1) { return 4;  }
    if constexpr (O == ElementOrder::O2) { return 10; }
    return 0;
}();

template<ElementOrder O> constexpr size_t EXTRA_NODES_IN_VOL_ELEM = [] {
    return NODES_IN_VOL_ELEM<O> - NODES_IN_VOL_ELEM<ElementOrder::O1>;
}();

template<> class Result<
    Phenomenon::ACOUSTIC,
    NumericalMethod::FEM,
    Domain::FREQUENCY,
    Dimension::D3
> {
public:
    Result();
    Result(const size_t&, const size_t&); // todo: make private
    ~Result();
    // todo: rule of 5
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> _data;
private:
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
    void set_frequency_range(
        const double&,
        const double&
    );
    void load_mesh(
        const char* const
    );
    void add_volume_material(
        const size_t&,
        const std::function<std::complex<double>(const double&)>&,
        const std::function<std::complex<double>(const double&)>&
    );
    void add_sound_source(
        const TypeOfSource&,
        const std::array<double,3>&,
        const PhysicalQuantity&,
        const std::function<std::complex<double>(const double&)>&
    );
    void add_sound_source(
        const TypeOfSource&,
        const size_t&,
        const PhysicalQuantity&,
        const std::function<std::complex<double>(const double&)>&
    );
    void add_surface_specific_acoustic_impedance(
        const size_t&,
        const std::function<std::complex<double>(const double&)>&
    );
    Result<
        Phenomenon::ACOUSTIC,
        NumericalMethod::FEM,
        Domain::FREQUENCY,
        Dimension::D3
    > run();

private:
    class Impl;
    std::unique_ptr<Impl> pimpl;
};

// alias for types
template<ElementOrder O>
using SimulationAcFemFreqD3 = typename numav::Simulation<
    Phenomenon::ACOUSTIC,
    NumericalMethod::FEM,
    Domain::FREQUENCY,
    Dimension::D3,
    O
>;
using ResultAcFemFreqD3 = typename numav::Result<
    Phenomenon::ACOUSTIC,
    NumericalMethod::FEM,
    Domain::FREQUENCY,
    Dimension::D3
>;

} // namespace numav
