// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include <functional>
#include <complex>
#include <memory>
#include <array>
#include <vector>

namespace numav {

enum class ElementOrder : size_t {
    O1,
    O2
};

enum class FrequencySamplingDensity : size_t {
    CONSTANT,
    QUADRATIC
};

#define NUMAV_ADD_VOL_MAT(T1, T2) \
    void add_volume_material( \
        const size_t external_volume_physical_group, \
        T1 density, \
        T2 sound_speed \
    );

#define NUMAV_ADD_SFC_MATERIAL(T) \
    void add_surface_material( \
        const size_t external_surface_physical_group, \
        const PhysicalQuantity physical_quantity_type, \
        T physical_quantity_value \
    );

#define NUMAV_ADD_SOUND_SOURCE(T1, T2) \
    void add_sound_source( \
        const TypeOfSource source_type, \
        T1 position, \
        const PhysicalQuantity physical_quantity_type, \
        T2 physical_quantity_value \
    );

#define NUMAV_ARRAY3 const std::array<Float,3UL>

#define NUMAV_DECLARE_SIM_AC_FEM_FREQ_D3_PUBLIC_METHODS \
    void set_maximum_frequency( \
        const Float maximum_frequnecy \
    ); \
    void set_frequency_range( \
        const Float minimum_frequnecy, \
        const Float maximum_frequnecy \
    ); \
    void set_frequency_steps_count( \
        const size_t frequency_steps_count \
    ); \
    void set_frequency_sampling_density( \
        const FrequencySamplingDensity frequency_sampling_density \
    ); \
    void set_frequency_steps( \
        const std::vector<Float>& frequency_steps \
    ); \
    void load_mesh( \
        const char* const path_to_mesh_file \
    ); \
    NUMAV_ADD_VOL_MAT(const FuncFloatToCmplx&, const FuncFloatToCmplx&)\
    NUMAV_ADD_VOL_MAT(const FuncFloatToCmplx&, const Cmplx            )\
    NUMAV_ADD_VOL_MAT(const FuncFloatToCmplx&, const char* const      )\
    NUMAV_ADD_VOL_MAT(const Cmplx            , const FuncFloatToCmplx&)\
    NUMAV_ADD_VOL_MAT(const Cmplx            , const Cmplx            )\
    NUMAV_ADD_VOL_MAT(const Cmplx            , const char* const      )\
    NUMAV_ADD_VOL_MAT(const char* const      , const FuncFloatToCmplx&)\
    NUMAV_ADD_VOL_MAT(const char* const      , const Cmplx            )\
    NUMAV_ADD_VOL_MAT(const char* const      , const char* const      )\
    NUMAV_ADD_SFC_MATERIAL(const FuncFloatToCmplx&) \
    NUMAV_ADD_SFC_MATERIAL(const Cmplx            ) \
    NUMAV_ADD_SFC_MATERIAL(const char* const      ) \
    NUMAV_ADD_SOUND_SOURCE(NUMAV_ARRAY3, const FuncFloatToCmplx&) \
    NUMAV_ADD_SOUND_SOURCE(NUMAV_ARRAY3, const Cmplx            ) \
    NUMAV_ADD_SOUND_SOURCE(NUMAV_ARRAY3, const char* const      ) \
    NUMAV_ADD_SOUND_SOURCE(const size_t, const FuncFloatToCmplx&) \
    NUMAV_ADD_SOUND_SOURCE(const size_t, const Cmplx            ) \
    NUMAV_ADD_SOUND_SOURCE(const size_t, const char* const      ) \
    void add_receiver( \
        const std::array<Float,3UL> coordinates \
    ); \
    void set_result_export_path( \
        const char* const path_to_nmvr_file \
    ); \
    void run( \
    );

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
    Simulation(const Simulation&) = delete;
    Simulation& operator=(const Simulation&) = delete;
    Simulation(Simulation&&) noexcept;
    Simulation& operator=(Simulation&&) noexcept;
    
    NUMAV_DECLARE_SIM_AC_FEM_FREQ_D3_PUBLIC_METHODS

private:
    class Impl;
    std::unique_ptr<Impl> _pimpl;
};

// alias for simulation type
template<ElementOrder O>
using SimulationAcFemFreqD3 = Simulation<
    Phenomenon::ACOUSTIC,
    NumericalMethod::FEM,
    Domain::FREQUENCY,
    Dimension::D3,
    O
>;

} // namespace numav
