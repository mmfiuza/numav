// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav/numav.hpp"
#include "modules/ac-fem-freq-d3/simulation/impl/impl.hpp"

namespace numav {

template <ElementOrder O>
SimulationAcFemFreqD3<O>::Simulation() {
    pimpl = std::make_unique<Impl>();
}

template <ElementOrder O>
SimulationAcFemFreqD3<O>::~Simulation() = default;

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