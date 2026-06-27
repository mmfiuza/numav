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
}

// Dimension enum class
JLCXX_MODULE define_module_Dimension(jlcxx::Module& mod) {
    mod.add_bits<Dimension>("type", jlcxx::julia_type("CppEnum"));
    mod.set_const("d3", Dimension::D3);
}

// SourceType enum class
JLCXX_MODULE define_module_SourceType(jlcxx::Module& mod) {
    mod.add_bits<SourceType>("type", jlcxx::julia_type("CppEnum"));
    mod.set_const("point", SourceType::POINT);
    mod.set_const("surface", SourceType::SURFACE);
}

// PhysicalQuantity enum class
JLCXX_MODULE define_module_PhysicalQuantity(jlcxx::Module& mod) {
    mod.add_bits<PhysicalQuantity>("type", jlcxx::julia_type("CppEnum"));
    mod.set_const("volume_velocity", PhysicalQuantity::VOLUME_VELOCITY);
    mod.set_const("particle_velocity", PhysicalQuantity::PARTICLE_VELOCITY);
    mod.set_const("pressure", PhysicalQuantity::PRESSURE);
    mod.set_const("impedance", PhysicalQuantity::IMPEDANCE);
}

// ElementShape enum class
JLCXX_MODULE define_module_ElementShape(jlcxx::Module& mod) {
    mod.add_bits<ElementShape>("type", jlcxx::julia_type("CppEnum"));
    mod.set_const("tetrahedron", ElementShape::TETRAHEDRON);
}

// ElementOrder enum class
JLCXX_MODULE define_module_ElementOrder(jlcxx::Module& mod) {
    mod.add_bits<ElementOrder>("type", jlcxx::julia_type("CppEnum"));
    mod.set_const("o1", ElementOrder::O1);
    mod.set_const("o2", ElementOrder::O2);
}

// FrequencySamplingDensity enum class
JLCXX_MODULE define_module_FrequencySamplingDensity(jlcxx::Module& mod) {
    mod.add_bits<FrequencySamplingDensity>("type",jlcxx::julia_type("CppEnum"));
    mod.set_const("constant", FrequencySamplingDensity::CONSTANT);
    mod.set_const("quadratic", FrequencySamplingDensity::QUADRATIC);
}

namespace jlcxx
{
    template<
        Phenomenon PHE,
        NumericalMethod NUM,
        Domain DOM,
        Dimension DIM,
        ElementShape SHA,
        ElementOrder ORD
    >
    struct BuildParameterList<Simulation<PHE, NUM, DOM, DIM, SHA, ORD>> {
        typedef ParameterList<
            std::integral_constant<Phenomenon, PHE>,
            std::integral_constant<NumericalMethod, NUM>,
            std::integral_constant<Domain, DOM>,
            std::integral_constant<Dimension, DIM>,
            std::integral_constant<ElementShape, SHA>,
            std::integral_constant<ElementOrder, ORD>
        > type;
    };
}

FuncFloatToCmplx cmplx_func(
    const jlcxx::SafeCFunction& real,
    const jlcxx::SafeCFunction& imag
) {
    auto f_real = jlcxx::make_function_pointer<Float(Float)>(real);
    auto f_imag = jlcxx::make_function_pointer<Float(Float)>(imag);
    return [f_real, f_imag](const Float f) -> Cmplx {
        return Cmplx(f_real(f), f_imag(f));
    };
}

std::array<Float,3UL> jl2cpp_array(const jlcxx::ArrayRef<Float> jl_array) {
    return std::array<Float,3UL>{jl_array[0], jl_array[1], jl_array[2]};
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    mod.add_type<
        jlcxx::Parametric<
            jlcxx::TypeVar<1>,
            jlcxx::TypeVar<2>,
            jlcxx::TypeVar<3>,
            jlcxx::TypeVar<4>,
            jlcxx::TypeVar<5>,
            jlcxx::TypeVar<6>
        >
    >("Simulation").apply<
        Simulation<
            Phenomenon::ACOUSTIC,
            NumericalMethod::FEM,
            Domain::FREQUENCY,
            Dimension::D3,
            ElementShape::TETRAHEDRON,
            ElementOrder::O1
        >,
        Simulation<
            Phenomenon::ACOUSTIC,
            NumericalMethod::FEM,
            Domain::FREQUENCY,
            Dimension::D3,
            ElementShape::TETRAHEDRON,
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
                const uint64_t freq_steps_count
            ) { w.set_frequency_steps_count(freq_steps_count); }
        );
        wrapped.module().method("set_frequency_sampling_density",
            []( WrappedT& w,
                const FrequencySamplingDensity pqv1
            ) { w.set_frequency_sampling_density(pqv1); }
        );
        wrapped.module().method("set_frequency_steps",
            []( WrappedT& w,
                const jlcxx::ArrayRef<Float> freq_steps
            ) { 
                std::vector<Float> cpp_vec;
                cpp_vec.reserve(freq_steps.size());
                for (auto& f : freq_steps) {
                    cpp_vec.emplace_back(f);
                }
                w.set_frequency_steps(cpp_vec);
            }
        );
        wrapped.module().method("load_mesh",
            []( WrappedT& w,
                const char* const path_to_mesh
            ) { w.load_mesh(path_to_mesh); }
        );

        // _add_volume_material
        wrapped.module().method("_add_volume_material",
            []( WrappedT& w,
                const uint64_t evpg,
                jlcxx::SafeCFunction pqv1_re,
                jlcxx::SafeCFunction pqv1_im,
                jlcxx::SafeCFunction pqv2_re,
                jlcxx::SafeCFunction pqv2_im
            ) { 
                w.add_volume_material(
                    evpg,
                    cmplx_func(pqv1_re, pqv1_im),
                    cmplx_func(pqv2_re, pqv2_im)
                );
            }
        );
        wrapped.module().method("_add_volume_material",
            []( WrappedT& w,
                const uint64_t evpg,
                jlcxx::SafeCFunction pqv1_re,
                jlcxx::SafeCFunction pqv1_im,
                const Cmplx pqv2
            ) { 
                w.add_volume_material(
                    evpg, cmplx_func(pqv1_re, pqv1_im), pqv2
                );
            }
        );
        wrapped.module().method("_add_volume_material",
            []( WrappedT& w,
                const uint64_t evpg,
                jlcxx::SafeCFunction pqv1_re,
                jlcxx::SafeCFunction pqv1_im,
                const char* const pqv2
            ) { 
                w.add_volume_material(
                    evpg, cmplx_func(pqv1_re, pqv1_im), pqv2
                );
            }
        );
        wrapped.module().method("_add_volume_material",
            []( WrappedT& w,
                const uint64_t evpg,
                const Cmplx pqv1,
                jlcxx::SafeCFunction pqv2_re,
                jlcxx::SafeCFunction pqv2_im
            ) { 
                w.add_volume_material(
                    evpg, pqv1, cmplx_func(pqv2_re, pqv2_im)
                );
            }
        );
        wrapped.module().method("_add_volume_material",
            []( WrappedT& w,
                const uint64_t evpg,
                const Cmplx pqv1,
                const Cmplx pqv2
            ) { w.add_volume_material(evpg, pqv1, pqv2); }
        );
        wrapped.module().method("_add_volume_material",
            []( WrappedT& w,
                const uint64_t evpg,
                const Cmplx pqv1,
                const char* const pqv2
            ) { w.add_volume_material(evpg, pqv1, pqv2); }
        );
        wrapped.module().method("_add_volume_material",
            []( WrappedT& w,
                const uint64_t evpg,
                const char* const pqv1,
                jlcxx::SafeCFunction pqv2_re,
                jlcxx::SafeCFunction pqv2_im
            ) { 
                w.add_volume_material(
                    evpg, pqv1, cmplx_func(pqv2_re, pqv2_im));
                }
        );
        wrapped.module().method("_add_volume_material",
            []( WrappedT& w,
                const uint64_t evpg,
                const char* const pqv1,
                const Cmplx pqv2
            ) { w.add_volume_material(evpg, pqv1, pqv2); }
        );
        wrapped.module().method("_add_volume_material",
            []( WrappedT& w,
                const uint64_t evpg,
                const char* const pqv1,
                const char* const pqv2
            ) { w.add_volume_material(evpg, pqv1, pqv2); }
        );

        // _add_surface_material
        wrapped.module().method("_add_surface_material",
            []( WrappedT& w,
                const uint64_t espg,
                const PhysicalQuantity pq,
                jlcxx::SafeCFunction pqv_re,
                jlcxx::SafeCFunction pqv_im
            ) { w.add_surface_material(espg, pq, cmplx_func(pqv_re, pqv_im)); }
        );
        wrapped.module().method("_add_surface_material",
            []( WrappedT& w,
                const uint64_t espg,
                const PhysicalQuantity pq,
                const Cmplx pqv
            ) { w.add_surface_material(espg, pq, pqv); }
        );
        wrapped.module().method("_add_surface_material",
            []( WrappedT& w,
                const uint64_t espg,
                const PhysicalQuantity pq,
                const char* const pqv
            ) { w.add_surface_material(espg, pq, pqv); }
        );
        
        // _add_sound_source
        wrapped.module().method("_add_sound_source",
            []( WrappedT& w,
                const SourceType source_type,
                const jlcxx::ArrayRef<Float> coords,
                const PhysicalQuantity pq,
                jlcxx::SafeCFunction pqv_re,
                jlcxx::SafeCFunction pqv_im
            ) { 
                w.add_sound_source(
                    source_type,
                    jl2cpp_array(coords),
                    pq, 
                    cmplx_func(pqv_re, pqv_im)
                );
            }
        );
        wrapped.module().method("_add_sound_source",
            []( WrappedT& w,
                const SourceType source_type,
                const jlcxx::ArrayRef<Float> coords,
                const PhysicalQuantity pq,
                const Cmplx pqv
            ) { 
                w.add_sound_source(source_type, jl2cpp_array(coords), pq, pqv);
            }
        );
        wrapped.module().method("_add_sound_source",
            []( WrappedT& w,
                const SourceType source_type,
                const jlcxx::ArrayRef<Float> coords,
                const PhysicalQuantity pq,
                const char* const pqv
            ) { 
                w.add_sound_source(source_type, jl2cpp_array(coords), pq, pqv);
            }
        );
        wrapped.module().method("_add_sound_source",
            []( WrappedT& w,
                const SourceType source_type,
                const uint64_t espg,
                const PhysicalQuantity pq,
                jlcxx::SafeCFunction pqv_re,
                jlcxx::SafeCFunction pqv_im
            ) { 
                w.add_sound_source(
                    source_type, espg, pq, cmplx_func(pqv_re, pqv_im)
                );
            }
        );
        wrapped.module().method("_add_sound_source",
            []( WrappedT& w,
                const SourceType source_type,
                const uint64_t espg,
                const PhysicalQuantity pq,
                const Cmplx pqv
            ) { w.add_sound_source(source_type, espg, pq, pqv); }
        );
        wrapped.module().method("_add_sound_source",
            []( WrappedT& w,
                const SourceType source_type,
                const uint64_t espg,
                const PhysicalQuantity pq,
                const char* const pqv
            ) { w.add_sound_source(source_type, espg, pq, pqv); }
        );

        wrapped.module().method("add_receiver", 
            []( WrappedT& w,
                const jlcxx::ArrayRef<Float> coords
            ) { w.add_receiver(jl2cpp_array(coords)); }
        );
        wrapped.module().method("set_result_export_path", 
            []( WrappedT& w,
                const char* const path_to_hdf5_file
            ) { w.set_result_export_path(path_to_hdf5_file); }
        );
        // "run" is a native Julia function, so we choose simulate as a name
        wrapped.module().method("simulate", 
            [](WrappedT& w) { w.run(); }
        );
    });
}
