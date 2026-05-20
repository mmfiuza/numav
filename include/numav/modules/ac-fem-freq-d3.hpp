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
    
    void set_maximum_frequency(
        const Float maximum_frequnecy
    );
    void set_frequency_range(
        const Float minimum_frequnecy,
        const Float maximum_frequnecy
    );
    void set_frequency_steps_count(
        const size_t frequency_steps_count
    );
    void set_frequency_sampling_density(
        const FrequencySamplingDensity frequency_sampling_density
    );
    void set_frequency_steps(
        const std::vector<Float>& frequency_steps
    );
    void load_mesh(
        const char* const path_to_mesh_file
    );
    void add_volume_material(
        const size_t external_volume_physical_group,
        const FuncFloatToCmplx& density,
        const FuncFloatToCmplx& sound_speed
    );
    void add_volume_material(
        const size_t external_volume_physical_group,
        const char* const density,
        const FuncFloatToCmplx& sound_speed
    );
    void add_volume_material(
        const size_t external_volume_physical_group,
        const FuncFloatToCmplx& density,
        const char* const sound_speed
    );
    void add_volume_material(
        const size_t external_volume_physical_group,
        const char* const density,
        const char* const sound_speed
    );
    void add_sound_source(
        const TypeOfSource source_type,
        const std::array<Float,3UL> coordinates,
        const PhysicalQuantity physical_quantity_type,
        const FuncFloatToCmplx& physical_quantity_value
    );
    void add_sound_source(
        const TypeOfSource source_type,
        const std::array<Float,3UL> coordinates,
        const PhysicalQuantity physical_quantity_type,
        const char* const physical_quantity_value
    );
    void add_sound_source(
        const TypeOfSource source_type,
        const size_t external_surface_physical_group,
        const PhysicalQuantity physical_quantity_type,
        const FuncFloatToCmplx& physical_quantity_value
    );
    void add_sound_source(
        const TypeOfSource source_type,
        const size_t external_surface_physical_group,
        const PhysicalQuantity physical_quantity_type,
        const char* const physical_quantity_value
    );
    void add_receiver(
        const std::array<Float,3UL> coordinates
    );
    void add_surface_material(
        const size_t external_surface_physical_group,
        const PhysicalQuantity physical_quantity_type,
        const FuncFloatToCmplx& physical_quantity_value
    );
    void add_surface_material(
        const size_t external_surface_physical_group,
        const PhysicalQuantity physical_quantity_type,
        const char* const physical_quantity_value
    );
    void set_result_export_path(
        const char* const path_to_nmvr_file
    );
    void run(
    );

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
