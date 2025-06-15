// Copyright (c) 2025 Matheus Machado Fiuza <matheus.fiuza@eac.ufsm.br>

#include "numav.hpp"
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"

// Phenomenon enum class
JLCXX_MODULE define_module_Phenomenon(jlcxx::Module& mod) {
    mod.add_bits<numav::Phenomenon>("Phenomenon_type", jlcxx::julia_type("CppEnum"));
    mod.set_const("acoustic", numav::Phenomenon::acoustic);
}

// NumericalMethod enum class
JLCXX_MODULE define_module_NumericalMethod(jlcxx::Module& mod) {
    mod.add_bits<numav::NumericalMethod>("NumericalMethod_type", jlcxx::julia_type("CppEnum"));
    mod.set_const("fem", numav::NumericalMethod::fem);
}

// Domain enum class
JLCXX_MODULE define_module_Domain(jlcxx::Module& mod) {
    mod.add_bits<numav::Domain>("Domain_type", jlcxx::julia_type("CppEnum"));
    mod.set_const("frequency", numav::Domain::frequency);
    mod.set_const("time", numav::Domain::time);
}

// Dimension enum class
JLCXX_MODULE define_module_Dimension(jlcxx::Module& mod) {
    mod.add_bits<numav::Dimension>("Dimension_type", jlcxx::julia_type("CppEnum"));
    mod.set_const("d1", numav::Dimension::d1);
    mod.set_const("d2", numav::Dimension::d2);
    mod.set_const("d3", numav::Dimension::d3);
}

// JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
// {   
//     mod.add_type<
//         jlcxx::Parametric<
//             jlcxx::TypeVar<1>, jlcxx::TypeVar<2>, jlcxx::TypeVar<3>, jlcxx::TypeVar<4>
//         >
//     >("Simulation").apply<
//         numav::Simulation<
//             numav::Phenomenon::acoustic,
//             numav::NumericalMethod::fem,
//             numav::Domain::frequency,
//             numav::Dimension::d3
//         >
//     > (
//         [](auto wrapped) {
//             typedef typename decltype(wrapped)::type WrappedT;
//             wrapped.template constructor<>();
//         }
//     );
// }



JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{   
    // mod.add_type<
    //     numav::Simulation<
    //         numav::Phenomenon::acoustic,
    //         numav::NumericalMethod::fem,
    //         numav::Domain::frequency,
    //         numav::Dimension::d3
    //     >
    // >("Simulation")
    // .constructor<>();


    // Helper to wrap NonTypeParam instances
    struct WrapNonTypeParam
    {
        template<typename TypeWrapperT>
        void operator()(TypeWrapperT&& wrapped)
        {
            typedef typename TypeWrapperT::type WrappedT;
            wrapped.template constructor<typename WrappedT::type>();
            // Access the module to add a free function
            // wrapped.module().method("get_nontype", [](const WrappedT& w) { return w.i; });
        }
    };

    mod.add_type<
        jlcxx::Parametric<
            jlcxx::TypeVar<1>, jlcxx::TypeVar<2>, jlcxx::TypeVar<3>, jlcxx::TypeVar<4>
        >
    >("Simulation").apply<
        numav::Simulation<
            numav::Phenomenon::acoustic,
            numav::NumericalMethod::fem,
            numav::Domain::frequency,
            numav::Dimension::d3
        >
    > (
        // [](auto wrapped) {
        //     typedef typename jlcxx::TypeWrapperT::type WrappedT;
        //     wrapped.template constructor<typename WrappedT::type>();
        // }
        WrapNonTypeParam()
    );
}