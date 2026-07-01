// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav/numav.hpp"
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"
#include "jlcxx/functions.hpp"

using namespace numav;

// NumericalMethod enum class
JLCXX_MODULE NumericalMethod_module(jlcxx::Module& mod) {
    mod.add_enum<NumericalMethod>("NumericalMethod",
        std::vector<const char*>({
            "fem"
        }),
        std::vector<uint64_t>({
            static_cast<uint64_t>(NumericalMethod::FEM)
        })
    );
}

// Equation enum class
JLCXX_MODULE Equation_module(jlcxx::Module& mod) {
    mod.add_enum<Equation>("Equation",
        std::vector<const char*>({
            "helmholtz"
        }),
        std::vector<uint64_t>({
            static_cast<uint64_t>(Equation::HELMHOLTZ)
        })
    );
}

// ElementShape enum class
JLCXX_MODULE ElementShape_module(jlcxx::Module& mod) {
    mod.add_enum<ElementShape>("ElementShape",
        std::vector<const char*>({
            "tetrahedron"
        }),
        std::vector<uint64_t>({
            static_cast<uint64_t>(ElementShape::TETRAHEDRON)
        })
    );
}

// ElementOrder enum class
JLCXX_MODULE ElementOrder_module(jlcxx::Module& mod) {
    mod.add_enum<ElementOrder>("ElementOrder",
        std::vector<const char*>({
            "linear",
            "quadratic"
        }),
        std::vector<uint64_t>({
            static_cast<uint64_t>(ElementOrder::LINEAR),
            static_cast<uint64_t>(ElementOrder::QUADRATIC)
        })
    );
}

// SourceType enum class
JLCXX_MODULE SourceType_module(jlcxx::Module& mod) {
    mod.add_enum<SourceType>("SourceType",
        std::vector<const char*>({
            "point",
            "surface"
        }),
        std::vector<uint64_t>({
            static_cast<uint64_t>(SourceType::POINT),
            static_cast<uint64_t>(SourceType::SURFACE),
        })
    );
}

// PhysicalQuantity enum class
JLCXX_MODULE PhysicalQuantity_module(jlcxx::Module& mod) {
    mod.add_enum<PhysicalQuantity>("PhysicalQuantity",
        std::vector<const char*>({
            "volume_velocity",
            "particle_velocity",
            "pressure",
            "impedance"
        }),
        std::vector<uint64_t>({
            static_cast<uint64_t>(PhysicalQuantity::VOLUME_VELOCITY),
            static_cast<uint64_t>(PhysicalQuantity::PARTICLE_VELOCITY),
            static_cast<uint64_t>(PhysicalQuantity::PRESSURE),
            static_cast<uint64_t>(PhysicalQuantity::IMPEDANCE)
        })
    );
}

// FrequencySamplingDensity enum class
JLCXX_MODULE FrequencySamplingDensity_module(jlcxx::Module& mod) {
    mod.add_enum<FrequencySamplingDensity>("FrequencySamplingDensity",
        std::vector<const char*>({
            "constant",
            "quadratic"
        }),
        std::vector<uint64_t>({
            static_cast<uint64_t>(FrequencySamplingDensity::CONSTANT),
            static_cast<uint64_t>(FrequencySamplingDensity::QUADRATIC)
        })
    );
}

namespace jlcxx
{
    template<NumericalMethod N, Equation E, ElementShape S, ElementOrder O>
    struct BuildParameterList<Simulation<N, E, S, O>> {
        typedef ParameterList<
            std::integral_constant<NumericalMethod, N>,
            std::integral_constant<Equation, E>,
            std::integral_constant<ElementShape, S>,
            std::integral_constant<ElementOrder, O>
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
            jlcxx::TypeVar<4>
        >
    >("Simulation").apply<
        Simulation<
            NumericalMethod::FEM,
            Equation::HELMHOLTZ,
            ElementShape::TETRAHEDRON,
            ElementOrder::LINEAR
        >,
        Simulation<
            NumericalMethod::FEM,
            Equation::HELMHOLTZ,
            ElementShape::TETRAHEDRON,
            ElementOrder::QUADRATIC
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
        wrapped.module().method("_run", 
            [](WrappedT& w) { w.run(); }
        );
    });
}
