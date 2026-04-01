// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/impl.hpp"

#include <tuple>
#include <cmath>
#include <numbers>
#include <algorithm>
#include <limits>
#include <chrono>
#include <iomanip>

#include "common/exception.hpp"
#include "common/hash-functions.hpp"
#include "common/maths.hpp"
#include "common/utils.hpp"
#include "modules/ac-fem-freq-d3/onemkl-solver.hpp"
#include "modules/ac-fem-freq-d3/eigen-solver.hpp"

namespace numav {

template<ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_solve_systems()
{
    // allocate the result matrix
    _cmplx_pressure_amp = 
        Eigen::Matrix<_cmplx_t, Eigen::Dynamic, Eigen::Dynamic>(
            _ni_count, _freq_count
        );

    // start timer
    auto start_time = std::chrono::system_clock::now();
    auto time_t_start = std::chrono::system_clock::to_time_t(start_time);
    std::cout << "Solver started at: "
        << std::put_time(std::localtime(&time_t_start), "%H:%M:%S") << "\n";
    std::cout << "Solution progress: 0%\n";

    for (size_t fi=0; fi!=_freq_count; ++fi)
    {
        _a_vals.fill(_cmplx_t(0.0, 0.0));
        _b_vals.fill(_cmplx_t(0.0, 0.0));
        const double freq = _freq_steps[fi];
        const double omega = 2*std::numbers::pi*freq;
        const double omega_squared = std::pow(omega, 2);

        // add point volume velocity to b vector
        for (size_t pvni=0; pvni!=_pvni_count; ++pvni)
        {
            const _cmplx_t volvel = (_pvni_to_forc_fi_part[pvni])(freq);
            *_pvni_to_ptr_in_b[pvni] += _cmplx_t(0.0,-omega) * volvel;
        }

        // add surface velocity to b vector
        for (size_t ispgv=0; ispgv!=_ispgv_count; ++ispgv)
        {
            const _cmplx_t velocity = (_ispgv_to_velocity[ispgv])(freq);
            const _cmplx_t fd_part = _cmplx_t(0.0,-omega) * velocity;
            for (size_t fipi=0; fipi!=_ispgv_to_ptr_in_b[ispgv].size(); ++fipi)
            {
                *_ispgv_to_ptr_in_b[ispgv][fipi] +=
                    _ispgv_to_forc_fi_part[ispgv][fipi] * fd_part;
            }
        }

        // add damping matrix to a
        for (size_t ispgi=0; ispgi!=_ispgi_count; ++ispgi) {
            const _cmplx_t impedance_value = _ispgi_to_impedance[ispgi](freq);
            const _cmplx_t damp_fd_part = _cmplx_t(0.0,omega)/impedance_value;
            
            for (size_t fipi=0; fipi!=_ispgi_to_ptr_in_a[ispgi].size(); ++fipi)
            {
                *_ispgi_to_ptr_in_a[ispgi][fipi] +=
                    _ispgi_to_damp_fi_part[ispgi][fipi] * damp_fd_part;
            }
        }

        // add stiffness and mass matrix to a
        for (size_t ivpg=0; ivpg!=_ivpg_count; ++ivpg) {
            const _cmplx_t density_value =
                (_ivpg_to_volprop[ivpg].density)(freq);
            const _cmplx_t soundspeed_value =
                (_ivpg_to_volprop[ivpg].soundspeed)(freq);

            const _cmplx_t stif_fd_part = 1.0 / density_value;
            const _cmplx_t mass_fd_part =
                - omega_squared / (density_value*std::pow(soundspeed_value,2));
            
            for (size_t fipi=0; fipi!=_ivpg_to_ptr_in_a[ivpg].size(); ++fipi) {
                *_ivpg_to_ptr_in_a[ivpg][fipi] +=
                    _ivpg_to_stif_fi_part[ivpg][fipi] * stif_fd_part +
                    _ivpg_to_mass_fi_part[ivpg][fipi] * mass_fd_part;
            }
        }

        // add pressure to a and b
        for (size_t pvi=0; pvi!=_pvi_count; ++pvi) {
            const _cmplx_t pressure = (_pvi_to_pressure[pvi])(freq);
            for (size_t fipi=0; fipi!=_pvi_to_ptr_in_a[pvi].size(); ++fipi) {
                *_pvi_to_ptr_in_a[pvi][fipi] += PENALTY_METHOD_CONSTANT;
                *_pvi_to_ptr_in_b[pvi][fipi] +=
                    PENALTY_METHOD_CONSTANT * pressure;
            }
        }

        // solve
        #if NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
            solve_using_eigen(
                _a_vals, _nnz_rowcol_idx_pairs, _b_vals, _b_row_idx,
                _cmplx_pressure_amp._data.data() + fi*_ni_count
            );
        #elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
            solve_using_onemkl(
                _dss_handle, _a_vals, _b_vals, _b_row_idx, _b_dense,
                _cmplx_pressure_amp.data() + fi*_ni_count
            );
        #else
            static_assert(false, "Invalid NUMAV_SYSTEM_SOLVER.");
        #endif

        // print the progress
        std::cout << "\033[A"; // Moves the cursor up one line
        std::cout << "\033[2K"; // Clears the entire line
        double progress_percentage =
            100.0*static_cast<double>(fi+1)/static_cast<double>(_freq_count);
        std::cout << "Solution progress: " << std::fixed << 
            std::setprecision(2) << progress_percentage << "%";

        if (fi > 99) { // print the time left // todo: fix the time
            auto current_time = std::chrono::system_clock::now();
            auto time_taken_until_now =
                std::chrono::duration_cast<std::chrono::seconds>(
                    current_time - start_time
                );
            auto time_left = 
                time_taken_until_now.count()/(fi+1)*(_freq_count-(fi+1));
            
            auto hours = time_left / 3600;
            auto minutes = (time_left % 3600) / 60;
            std::cout << " - estimated time left: " <<
                hours << " h and " << minutes << " min\n";
        }
        else {
            std::cout << "\n";
        }
    }
    _did_run = true;
    #if NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
        constexpr MKL_INT options = NUMAV_MKL_OPTIONS;
        _INTEGER_t error_id = dss_delete(_dss_handle, options);
        if (error_id != MKL_DSS_SUCCESS) { print_dss_error(error_id); }
    #endif

    // print finish
    auto end_time = std::chrono::system_clock::now();
    auto time_t_end = std::chrono::system_clock::to_time_t(end_time);
    std::cout << "Solver ended at: " << std::put_time(std::localtime(&time_t_end), "%H:%M:%S") << "\n";
}

// explicit instantiation declarations
INSTANTIATE_SIMULATION_CLASS

} // namespace numav
