// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

namespace numav {

enum class Equation : uint64_t {
    HELMHOLTZ
};

enum class ElementShape : uint64_t {
    TETRAHEDRON
};

enum class ElementOrder : uint64_t {
    LINEAR,
    QUADRATIC
};

enum class SourceType : uint64_t {
    POINT,
    SURFACE
};

enum class PhysicalQuantity : uint64_t {
    VOLUME_VELOCITY,
    PARTICLE_VELOCITY,
    PRESSURE,
    IMPEDANCE
};

enum class FrequencySamplingDensity : uint64_t {
    CONSTANT,
    QUADRATIC
};

#define NUMAV_PUBLIC_METHOD(method_name, ...) void method_name(__VA_ARGS__);

#define NUMAV_SIM_AC_FEM_FREQ_D3_PUBLIC_METHODS\
    NUMAV_PUBLIC_METHOD(set_maximum_frequency,\
        const Float\
    )\
    NUMAV_PUBLIC_METHOD(set_frequency_range,\
        const Float, const Float\
    )\
    NUMAV_PUBLIC_METHOD(set_frequency_steps_count,\
        const uint64_t\
    )\
    NUMAV_PUBLIC_METHOD(set_frequency_sampling_density,\
        const FrequencySamplingDensity\
    )\
    NUMAV_PUBLIC_METHOD(set_frequency_steps,\
        const std::vector<Float>&\
    )\
    NUMAV_PUBLIC_METHOD(load_mesh,\
        const char* const\
    )\
    NUMAV_PUBLIC_METHOD(add_volume_material,\
        const uint64_t, const FuncFloatToCmplx&, const FuncFloatToCmplx&\
    )\
    NUMAV_PUBLIC_METHOD(add_volume_material,\
        const uint64_t, const FuncFloatToCmplx&, const Cmplx\
    )\
    NUMAV_PUBLIC_METHOD(add_volume_material,\
        const uint64_t, const FuncFloatToCmplx&, const char* const\
    )\
    NUMAV_PUBLIC_METHOD(add_volume_material,\
        const uint64_t, const Cmplx, const FuncFloatToCmplx&\
    )\
    NUMAV_PUBLIC_METHOD(add_volume_material,\
        const uint64_t, const Cmplx, const Cmplx\
    )\
    NUMAV_PUBLIC_METHOD(add_volume_material,\
        const uint64_t, const Cmplx, const char* const\
    )\
    NUMAV_PUBLIC_METHOD(add_volume_material,\
        const uint64_t, const char* const, const FuncFloatToCmplx&\
    )\
    NUMAV_PUBLIC_METHOD(add_volume_material,\
        const uint64_t, const char* const, const Cmplx\
    )\
    NUMAV_PUBLIC_METHOD(add_volume_material,\
        const uint64_t, const char* const, const char* const\
    )\
    NUMAV_PUBLIC_METHOD(add_surface_material,\
        const uint64_t, const PhysicalQuantity, const FuncFloatToCmplx&\
    )\
    NUMAV_PUBLIC_METHOD(add_surface_material,\
        const uint64_t, const PhysicalQuantity, const Cmplx\
    )\
    NUMAV_PUBLIC_METHOD(add_surface_material,\
        const uint64_t, const PhysicalQuantity, const char* const\
    )\
    NUMAV_PUBLIC_METHOD(add_sound_source,\
        const SourceType, const Coord,\
        const PhysicalQuantity, const FuncFloatToCmplx&\
    )\
    NUMAV_PUBLIC_METHOD(add_sound_source,\
        const SourceType, const Coord,\
        const PhysicalQuantity, const Cmplx\
    )\
    NUMAV_PUBLIC_METHOD(add_sound_source,\
        const SourceType, const Coord,\
        const PhysicalQuantity, const char* const\
    )\
    NUMAV_PUBLIC_METHOD(add_sound_source,\
        const SourceType, const uint64_t,\
        const PhysicalQuantity, const FuncFloatToCmplx&\
    )\
    NUMAV_PUBLIC_METHOD(add_sound_source,\
        const SourceType, const uint64_t,\
        const PhysicalQuantity, const Cmplx\
    )\
    NUMAV_PUBLIC_METHOD(add_sound_source,\
        const SourceType, const uint64_t,\
        const PhysicalQuantity, const char* const\
    )\
    NUMAV_PUBLIC_METHOD(add_receiver,\
        const Coord\
    )\
    NUMAV_PUBLIC_METHOD(set_result_export_path,\
        const char* const\
    )\
    NUMAV_PUBLIC_METHOD(run\
    )

template<ElementOrder O>
class Simulation<
    NumericalMethod::FEM,
    Equation::HELMHOLTZ,
    ElementShape::TETRAHEDRON,
    O
> {
public:
    Simulation();
    ~Simulation();
    Simulation(const Simulation&) = delete;
    Simulation& operator=(const Simulation&) = delete;
    Simulation(Simulation&&) noexcept;
    Simulation& operator=(Simulation&&) noexcept;
    
    NUMAV_SIM_AC_FEM_FREQ_D3_PUBLIC_METHODS

private:
    class Impl;
    std::unique_ptr<Impl> _pimpl;
};

// alias for simulation type
template<ElementOrder O>
using SimulationFemHelmTet = Simulation<
    NumericalMethod::FEM,
    Equation::HELMHOLTZ,
    ElementShape::TETRAHEDRON,
    O
>;

} // namespace numav
