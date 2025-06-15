// Copyright (c) 2025 Matheus Machado Fiuza <matheus.fiuza@eac.ufsm.br>

#include "numav.hpp"
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"

// Phenomenon enum class
JLCXX_MODULE define_module_Phenomenon(jlcxx::Module& mod) {
    mod.add_bits<numav::Phenomenon>("Phenomenon_type", jlcxx::julia_type("CppEnum"));
    mod.set_const("acoustic", numav::Phenomenon::ACOUSTIC);
}

// NumericalMethod enum class
JLCXX_MODULE define_module_NumericalMethod(jlcxx::Module& mod) {
    mod.add_bits<numav::NumericalMethod>("NumericalMethod_type", jlcxx::julia_type("CppEnum"));
    mod.set_const("fem", numav::NumericalMethod::FEM);
}

// Domain enum class
JLCXX_MODULE define_module_Domain(jlcxx::Module& mod) {
    mod.add_bits<numav::Domain>("Domain_type", jlcxx::julia_type("CppEnum"));
    mod.set_const("frequency", numav::Domain::FREQUENCY);
    mod.set_const("time", numav::Domain::TIME);
}

// Dimension enum class
JLCXX_MODULE define_module_Dimension(jlcxx::Module& mod) {
    mod.add_bits<numav::Dimension>("Dimension_type", jlcxx::julia_type("CppEnum"));
    mod.set_const("d1", numav::Dimension::D1);
    mod.set_const("d2", numav::Dimension::D2);
    mod.set_const("d3", numav::Dimension::D3);
}

namespace jlcxx
{
    template<
    numav::Phenomenon PHENOMENON,
    numav::NumericalMethod NUMERICAL_METHOD,
    numav::Domain DOMAIN,
    numav::Dimension DIMENSION
    >
    struct BuildParameterList<numav::Simulation<
    PHENOMENON, NUMERICAL_METHOD, DOMAIN, DIMENSION
    >> {
        typedef ParameterList<
        std::integral_constant<numav::Phenomenon, PHENOMENON>,
        std::integral_constant<numav::NumericalMethod, NUMERICAL_METHOD>,
        std::integral_constant<numav::Domain, DOMAIN>,
        std::integral_constant<numav::Dimension, DIMENSION>
        > type;
    };
} // namespace jlcxx

struct WrapSimulation
{
    template<typename TypeWrapperT>
    void operator()(TypeWrapperT&& wrapped)
    {
        typedef typename TypeWrapperT::type WrappedT;
        wrapped.template constructor<>();
    }
};

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{   
    mod.add_type<
        jlcxx::Parametric<
            jlcxx::TypeVar<1>, jlcxx::TypeVar<2>, jlcxx::TypeVar<3>, jlcxx::TypeVar<4>
        >
    >("Simulation").apply<
        numav::Simulation<
            numav::Phenomenon::ACOUSTIC,
            numav::NumericalMethod::FEM,
            numav::Domain::FREQUENCY,
            numav::Dimension::D3
        >
    > (WrapSimulation());
}