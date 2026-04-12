// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include <functional>
#include <complex>
#include <memory>
#include <vector>

namespace numav {

enum class FrequencySamplingDensity {
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
        const double&
    );
    void set_frequency_range(
        const double&,
        const double&
    );
    void set_frequency_steps_count(
        const size_t&
    );
    void set_frequency_sampling_density(
        const FrequencySamplingDensity&
    );
    void set_frequency_steps(
        const std::vector<double>&
    );
    void load_mesh(
        const char* const
    );
    void add_volume_material(
        const size_t&,
        const std::function<std::complex<double>(const double&)>&,
        const std::function<std::complex<double>(const double&)>&
    );
    void add_volume_material(
        const size_t&,
        const char* const,
        const std::function<std::complex<double>(const double&)>&
    );
    void add_volume_material(
        const size_t&,
        const std::function<std::complex<double>(const double&)>&,
        const char* const
    );
    void add_volume_material(
        const size_t&,
        const char* const,
        const char* const
    );
    void add_sound_source(
        const TypeOfSource&,
        const std::array<double,3>&,
        const PhysicalQuantity&,
        const std::function<std::complex<double>(const double&)>&
    );
    void add_sound_source(
        const TypeOfSource&,
        const std::array<double,3>&,
        const PhysicalQuantity&,
        const char* const
    );
    void add_sound_source(
        const TypeOfSource&,
        const size_t&,
        const PhysicalQuantity&,
        const std::function<std::complex<double>(const double&)>&
    );
    void add_sound_source(
        const TypeOfSource&,
        const size_t&,
        const PhysicalQuantity&,
        const char* const
    );
    void add_surface_specific_acoustic_impedance(
        const size_t&,
        const std::function<std::complex<double>(const double&)>&
    );
    void add_surface_specific_acoustic_impedance(
        const size_t&,
        const char* const
    );
    void run(
    );
    void export_result(
        const char* const
    );

private:
    class Impl;
    std::unique_ptr<Impl> pimpl;
};

// alias for simulation type
template<ElementOrder O>
using SimulationAcFemFreqD3 = typename numav::Simulation<
    Phenomenon::ACOUSTIC,
    NumericalMethod::FEM,
    Domain::FREQUENCY,
    Dimension::D3,
    O
>;

} // namespace numav
