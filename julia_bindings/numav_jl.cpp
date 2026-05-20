// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav/numav.hpp"
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"
#include "jlcxx/functions.hpp"

using namespace numav;

// Phenomenon enum class
JLCXX_MODULE define_module_Phenomenon(jlcxx::Module& mod) {
    mod.add_bits<Phenomenon>("type", jlcxx::julia_type("CppEnum"));
    mod.set_const("acoustic", Phenomenon::ACOUSTIC);
}

// NumericalMethod enum class
JLCXX_MODULE define_module_NumericalMethod(jlcxx::Module& mod) {
    mod.add_bits<NumericalMethod>("type", jlcxx::julia_type("CppEnum"));
    mod.set_const("fem", NumericalMethod::FEM);
}

// Domain enum class
JLCXX_MODULE define_module_Domain(jlcxx::Module& mod) {
    mod.add_bits<Domain>("type", jlcxx::julia_type("CppEnum"));
    mod.set_const("frequency", Domain::FREQUENCY);
    mod.set_const("time", Domain::TIME);
}

// Dimension enum class
JLCXX_MODULE define_module_Dimension(jlcxx::Module& mod) {
    mod.add_bits<Dimension>("type", jlcxx::julia_type("CppEnum"));
    mod.set_const("d1", Dimension::D1);
    mod.set_const("d2", Dimension::D2);
    mod.set_const("d3", Dimension::D3);
}

// TypeOfSource enum class
JLCXX_MODULE define_module_TypeOfSource(jlcxx::Module& mod) {
    mod.add_bits<numav::TypeOfSource>("type", jlcxx::julia_type("CppEnum"));
    mod.set_const("point", numav::TypeOfSource::POINT);
    mod.set_const("surface", numav::TypeOfSource::SURFACE);
}

// PhysicalQuantity enum class
JLCXX_MODULE define_module_PhysicalQuantity(jlcxx::Module& mod) {
    mod.add_bits<numav::PhysicalQuantity>(
        "type", jlcxx::julia_type("CppEnum")
    );
    mod.set_const(
        "volume_velocity", numav::PhysicalQuantity::VOLUME_VELOCITY
    );
    mod.set_const(
        "particle_velocity", numav::PhysicalQuantity::PARTICLE_VELOCITY
    );
    mod.set_const(
        "pressure", numav::PhysicalQuantity::PRESSURE
    );
    mod.set_const(
        "impedance", numav::PhysicalQuantity::IMPEDANCE
    );
}

// ElementOrder enum class
JLCXX_MODULE define_module_ElementOrder(jlcxx::Module& mod) {
    mod.add_bits<ElementOrder>("type", jlcxx::julia_type("CppEnum"));
    mod.set_const("o1", ElementOrder::O1);
    mod.set_const("o2", ElementOrder::O2);
}

// FrequencySamplingDensity enum class
JLCXX_MODULE define_module_FrequencySamplingDensity(jlcxx::Module& mod) {
    mod.add_bits<numav::FrequencySamplingDensity>(
        "type", jlcxx::julia_type("CppEnum")
    );
    mod.set_const("constant", numav::FrequencySamplingDensity::CONSTANT);
    mod.set_const("quadratic", numav::FrequencySamplingDensity::QUADRATIC);
}

namespace jlcxx
{
    template<
        Phenomenon PHE,
        NumericalMethod NUM,
        Domain DOM,
        Dimension DIM,
        ElementOrder ORD
    >
    struct BuildParameterList<Simulation<PHE, NUM, DOM, DIM, ORD>> {
        typedef ParameterList<
            std::integral_constant<Phenomenon, PHE>,
            std::integral_constant<NumericalMethod, NUM>,
            std::integral_constant<Domain, DOM>,
            std::integral_constant<Dimension, DIM>,
            std::integral_constant<ElementOrder, ORD>
        > type;
    };
}

FuncFloatToCmplx cmplx_func(
    jlcxx::SafeCFunction real,
    jlcxx::SafeCFunction imag
) {
    auto f_real = jlcxx::make_function_pointer<Float(Float)>(real);
    auto f_imag = jlcxx::make_function_pointer<Float(Float)>(imag);
    return [f_real, f_imag](const Float f) -> Cmplx {
        return Cmplx(f_real(f), f_imag(f));
    };
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    mod.add_type<
        jlcxx::Parametric<
            jlcxx::TypeVar<1>,
            jlcxx::TypeVar<2>,
            jlcxx::TypeVar<3>,
            jlcxx::TypeVar<4>,
            jlcxx::TypeVar<5>
        >
    >("Simulation").apply<
        Simulation<
            Phenomenon::ACOUSTIC,
            NumericalMethod::FEM,
            Domain::FREQUENCY,
            Dimension::D3,
            ElementOrder::O1
        >,
        Simulation<
            Phenomenon::ACOUSTIC,
            NumericalMethod::FEM,
            Domain::FREQUENCY,
            Dimension::D3,
            ElementOrder::O2
        >
    >([](auto&& wrapped) {
        using WrappedT = typename std::decay_t<decltype(wrapped)>::type;
        wrapped.template constructor<>();
        wrapped.module().method("set_maximum_frequency",
            []( WrappedT& w,
                const Float freq_max
            ) { w.set_maximum_frequency(freq_max); }
        );
        wrapped.module().method("set_frequency_range",
            []( WrappedT& w,
                const Float freq_min,
                const Float freq_max
            ) { w.set_frequency_range(freq_min, freq_max); }
        );
        wrapped.module().method("set_frequency_steps_count",
            []( WrappedT& w,
                const size_t freq_steps_count
            ) { w.set_frequency_steps_count(freq_steps_count); }
        );
        wrapped.module().method("set_frequency_sampling_density",
            []( WrappedT& w,
                const numav::FrequencySamplingDensity density
            ) { w.set_frequency_sampling_density(density); }
        );
        wrapped.module().method("set_frequency_steps",
            []( WrappedT& w,
                const std::vector<Float>& freq_steps
            ) { w.set_frequency_steps(freq_steps); }
        );
        wrapped.module().method("load_mesh",
            []( WrappedT& w,
                const char* const path_to_mesh
            ) { w.load_mesh(path_to_mesh); }
        );
        wrapped.module().method("_add_volume_material",
            []( WrappedT& w,
                const size_t evpg,
                jlcxx::SafeCFunction density_func_real,
                jlcxx::SafeCFunction density_func_imag,
                jlcxx::SafeCFunction soundspeed_func_real,
                jlcxx::SafeCFunction soundspeed_func_imag
            ) { 
                w.add_volume_material(
                    evpg,
                    cmplx_func(density_func_real, density_func_imag),
                    cmplx_func(soundspeed_func_real, soundspeed_func_imag)
                );
            }
        );
        wrapped.module().method("_add_volume_material",
            []( WrappedT& w,
                const size_t evpg,
                const char* const rho_table,
                jlcxx::SafeCFunction soundspeed_func_real,
                jlcxx::SafeCFunction soundspeed_func_imag
            ) { 
                w.add_volume_material(
                    evpg,
                    rho_table,
                    cmplx_func(soundspeed_func_real, soundspeed_func_imag));
                }
        );
        wrapped.module().method("_add_volume_material",
            []( WrappedT& w,
                const size_t evpg,
                jlcxx::SafeCFunction density_func_real,
                jlcxx::SafeCFunction density_func_imag,
                const char* const c_table
            ) { 
                w.add_volume_material(
                    evpg,
                    cmplx_func(density_func_real, density_func_imag),
                    c_table
                );
            }
        );
        wrapped.module().method("_add_volume_material",
            []( WrappedT& w,
                const size_t evpg,
                const char* const rho_table,
                const char* const c_table
            ) { w.add_volume_material(evpg, rho_table, c_table); }
        );
        wrapped.module().method("_add_sound_source",
            []( WrappedT& w,
                const numav::TypeOfSource source_type,
                const jlcxx::ArrayRef<Float> coords,
                const numav::PhysicalQuantity pq_type,
                jlcxx::SafeCFunction pq_func_real,
                jlcxx::SafeCFunction pq_func_imag
            ) { 
                w.add_sound_source(
                    source_type,
                    std::array<Float,3UL>{coords[0], coords[1], coords[2]},
                    pq_type, 
                    cmplx_func(pq_func_real, pq_func_imag)
                );
            }
        );
        wrapped.module().method("_add_sound_source",
            []( WrappedT& w,
                const numav::TypeOfSource source_type,
                const jlcxx::ArrayRef<Float> coords,
                const numav::PhysicalQuantity pq_type,
                const char* const pq_table
            ) { 
                w.add_sound_source(
                    source_type,
                    std::array<Float,3UL>{coords[0], coords[1], coords[2]},
                    pq_type, 
                    pq_table
                );
            }
        );
        wrapped.module().method("_add_sound_source",
            []( WrappedT& w,
                const numav::TypeOfSource source_type,
                const size_t espg,
                const numav::PhysicalQuantity pq_type,
                jlcxx::SafeCFunction pq_func_real,
                jlcxx::SafeCFunction pq_func_imag
            ) { 
                w.add_sound_source(
                    source_type,
                    espg,
                    pq_type, 
                    cmplx_func(pq_func_real, pq_func_imag)
                );
            }
        );
        wrapped.module().method("_add_sound_source",
            []( WrappedT& w,
                const numav::TypeOfSource source_type,
                const size_t espg,
                const numav::PhysicalQuantity pq_type,
                const char* const pq_table
            ) { 
                w.add_sound_source(
                    source_type,
                    espg,
                    pq_type, 
                    pq_table
                );
            }
        );
        wrapped.module().method("_add_surface_material",
            []( WrappedT& w,
                const size_t espg,
                const PhysicalQuantity pq_type,
                jlcxx::SafeCFunction pq_func_real,
                jlcxx::SafeCFunction pq_func_imag
            ) {
                w.add_surface_material(
                    espg,
                    pq_type,
                    cmplx_func(pq_func_real, pq_func_imag)
                );
            }
        );
        wrapped.module().method("_add_surface_material",
            []( WrappedT& w,
                const size_t espg,
                const PhysicalQuantity pq_type,
                const char* const pq_table
            ) { 
                w.add_surface_material(
                    espg,
                    pq_type,
                    pq_table
                );
            }
        );
        wrapped.module().method("set_result_export_path", 
            []( WrappedT& w,
                const char* const path_to_nmvr_file
            ) { w.set_result_export_path(path_to_nmvr_file); }
        );
        // "run" is a native Julia function, so we choose simulate as a name
        wrapped.module().method("simulate", 
            [](WrappedT& w) { w.run(); }
        );
    });
}
