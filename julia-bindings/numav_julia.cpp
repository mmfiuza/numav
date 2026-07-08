// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav/numav.hpp"
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"
#include "jlcxx/functions.hpp"

using namespace numav;

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

std::array<Float,3UL> jl2cpp_array(const jlcxx::ArrayRef<Float> jl_array) {
    return std::array<Float,3UL>{jl_array[0], jl_array[1], jl_array[2]};
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    mod.add_enum<NumericalMethod>("_cpp_NumericalMethod",
        std::vector<const char*>({
            "_cpp_NumericalMethod_fem"
        }),
        std::vector<uint64_t>({
            static_cast<uint64_t>(NumericalMethod::FEM)
        })
    );
    mod.add_enum<Equation>("_cpp_Equation",
        std::vector<const char*>({
            "_cpp_Equation_helmholtz"
        }),
        std::vector<uint64_t>({
            static_cast<uint64_t>(Equation::HELMHOLTZ)
        })
    );
    mod.add_enum<ElementShape>("_cpp_ElementShape",
        std::vector<const char*>({
            "_cpp_ElementShape_tetrahedron"
        }),
        std::vector<uint64_t>({
            static_cast<uint64_t>(ElementShape::TETRAHEDRON)
        })
    );
    mod.add_enum<ElementOrder>("_cpp_ElementOrder",
        std::vector<const char*>({
            "_cpp_ElementOrder_linear",
            "_cpp_ElementOrder_quadratic"
        }),
        std::vector<uint64_t>({
            static_cast<uint64_t>(ElementOrder::LINEAR),
            static_cast<uint64_t>(ElementOrder::QUADRATIC)
        })
    );
    mod.add_enum<SourceType>("_cpp_SourceType",
        std::vector<const char*>({
            "_cpp_SourceType_point",
            "_cpp_SourceType_surface"
        }),
        std::vector<uint64_t>({
            static_cast<uint64_t>(SourceType::POINT),
            static_cast<uint64_t>(SourceType::SURFACE),
        })
    );
    mod.add_enum<PhysicalQuantity>("_cpp_PhysicalQuantity",
        std::vector<const char*>({
            "_cpp_PhysicalQuantity_volume_velocity",
            "_cpp_PhysicalQuantity_particle_velocity",
            "_cpp_PhysicalQuantity_pressure",
            "_cpp_PhysicalQuantity_impedance"
        }),
        std::vector<uint64_t>({
            static_cast<uint64_t>(PhysicalQuantity::VOLUME_VELOCITY),
            static_cast<uint64_t>(PhysicalQuantity::PARTICLE_VELOCITY),
            static_cast<uint64_t>(PhysicalQuantity::PRESSURE),
            static_cast<uint64_t>(PhysicalQuantity::IMPEDANCE)
        })
    );
    mod.add_type<
        jlcxx::Parametric<
            jlcxx::TypeVar<1>,
            jlcxx::TypeVar<2>,
            jlcxx::TypeVar<3>,
            jlcxx::TypeVar<4>
        >
    >("_cpp_Simulation").apply<
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
        wrapped.module().method("_cpp_set_frequency_vector!",
            []( WrappedT& w,
                const jlcxx::ArrayRef<Float> freq_steps
            ) { 
                std::vector<Float> cpp_vec;
                cpp_vec.reserve(freq_steps.size());
                for (auto& f : freq_steps) {
                    cpp_vec.emplace_back(f);
                }
                w.set_frequency_vector(cpp_vec);
            }
        );
        wrapped.module().method("_cpp_load_mesh!",
            []( WrappedT& w,
                const char* const path_to_mesh
            ) { w.load_mesh(path_to_mesh); }
        );
        wrapped.module().method("_cpp_add_volume_material!",
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
        wrapped.module().method("_cpp_add_surface_material!",
            []( WrappedT& w,
                const uint64_t espg,
                const PhysicalQuantity pq,
                jlcxx::SafeCFunction pqv_re,
                jlcxx::SafeCFunction pqv_im
            ) { w.add_surface_material(espg, pq, cmplx_func(pqv_re, pqv_im)); }
        );
        // _add_sound_source
        wrapped.module().method("_cpp_add_sound_source!",
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
        wrapped.module().method("_cpp_add_sound_source!",
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
        wrapped.module().method("_cpp_set_result_export_path!", 
            []( WrappedT& w,
                const char* const path_to_hdf5_file
            ) { w.set_result_export_path(path_to_hdf5_file); }
        );
        wrapped.module().method("_cpp_run!", 
            [](WrappedT& w) { w.run(); }
        );
    });
}
