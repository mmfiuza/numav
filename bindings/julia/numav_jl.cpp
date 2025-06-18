// Copyright (c) 2025 Matheus Machado Fiuza <matheus.fiuza@eac.ufsm.br>

#include "numav.hpp"
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"
#include "jlcxx/functions.hpp"
#include "julia.h"

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

// TypeOfSource enum class
JLCXX_MODULE define_module_TypeOfSource(jlcxx::Module& mod) {
    mod.add_bits<numav::TypeOfSource>("TypeOfSource_type", jlcxx::julia_type("CppEnum"));
    mod.set_const("point", numav::TypeOfSource::POINT);
    mod.set_const("surface", numav::TypeOfSource::SURFACE);
}

// PhysicalQuantity enum class
JLCXX_MODULE define_module_PhysicalQuantity(jlcxx::Module& mod) {
    mod.add_bits<numav::PhysicalQuantity>("PhysicalQuantity_type", jlcxx::julia_type("CppEnum"));
    mod.set_const("pressure", numav::PhysicalQuantity::PRESSURE);
    mod.set_const("volume_velocity", numav::PhysicalQuantity::VOLUME_VELOCITY);
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

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{   
    mod.add_type<jlcxx::Parametric<
        jlcxx::TypeVar<1>, jlcxx::TypeVar<2>, jlcxx::TypeVar<3>, jlcxx::TypeVar<4>
    >>("Simulation").apply<numav::Simulation<
        numav::Phenomenon::ACOUSTIC,
        numav::NumericalMethod::FEM,
        numav::Domain::FREQUENCY,
        numav::Dimension::D3
    >>(
        [](auto&& wrapped) {
            using WrappedT = typename std::decay_t<decltype(wrapped)>::type;
            wrapped.template constructor<>();
            wrapped.module().method("set_element_order",
                []( WrappedT& w,
                    const uint64_t& order
                ) { w.set_element_order(order); }
            );
            wrapped.module().method("set_freq_limits",
                []( WrappedT& w,
                    const double& freq_min,
                    const double& freq_max
                ) { w.set_freq_limits(freq_min, freq_max); }
            );
            wrapped.module().method("load_mesh",
                []( WrappedT& w,
                    const std::string& path_to_mesh
                ) { w.load_mesh(path_to_mesh); }
            );
            wrapped.module().method("add_volume_material",
                []( WrappedT& w,
                    const uint64_t& id,
                    const double& rho,
                    const double& c
                ) { w.add_volume_material(id, rho, c); }
            );
            wrapped.module().method("_add_source",
                []( WrappedT& w,
                    const numav::TypeOfSource& source_type,
                    const jlcxx::ArrayRef<double> p,
                    const numav::PhysicalQuantity& quantity_type,
                    jlcxx::SafeCFunction quantity_value
                ) { 
                    auto quantity_value_ptr = 
                        jlcxx::make_function_pointer<double(double)>(quantity_value);
                    w.add_source(source_type, std::array<double,3>{p[0],p[1],p[2]},
                        quantity_type, quantity_value_ptr);
                }
            );
        }
    );
}
