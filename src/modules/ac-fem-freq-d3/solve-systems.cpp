// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/impl.hpp"

#include <cmath>
#include <numbers>

#include "common/nmvr-format.hpp"
#include "modules/ac-fem-freq-d3/ldlt-solver-solver.hpp"
#include "modules/ac-fem-freq-d3/onemkl-solver.hpp"
#include "modules/ac-fem-freq-d3/eigen-solver.hpp"

#include <indicators/progress_bar.hpp>
#include <indicators/cursor_control.hpp>

namespace numav {

template<ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_solve_systems()
{
    _begin_nmvr_file();
    _write_simulation_inputs_to_nmvr_file();
    nmvr::write_data_chunk_header(
        _nvmr_file,
        nmvr::COMPLEX_AMPLITUDE_OF_SOUND_PRESSURE_SOLUTION_CHUNK_ID,
        _ni_count * _freq_count * sizeof(_cmplx_t)
    );

    // print start time
    auto start_time = std::chrono::system_clock::now();
    auto time_t_start = std::chrono::system_clock::to_time_t(start_time);
    std::cout << "Solver started at: "
        << std::put_time(std::localtime(&time_t_start), "%Hh:%Mm:%Ss") << "\n";

    // create progress bar
    indicators::show_console_cursor(false);
    indicators::ProgressBar bar {
        indicators::option::BarWidth{37UL},
        indicators::option::Start{" |"},
        indicators::option::Fill{"="},
        indicators::option::Lead{"="},
        indicators::option::Remainder{" "},
        indicators::option::End{"|"},
        indicators::option::PrefixText{"Running"},
        indicators::option::ForegroundColor{indicators::Color::unspecified},
        indicators::option::ShowPercentage{true},
        indicators::option::ShowElapsedTime{true},
        indicators::option::ShowRemainingTime{true},
        indicators::option::FontStyles{
            std::vector<indicators::FontStyle>{indicators::FontStyle::bold}
        },
        indicators::option::MinProgress{0UL},
        indicators::option::MaxProgress{_freq_count}
    };

    for (size_t fi = 0UL; fi != _freq_count; ++fi)
    {
        // reset values of A matrix and b vector
        _a_vals.fill(_cmplx_t(0.0, 0.0));
        _b_vals.fill(_cmplx_t(0.0, 0.0));
        #if NUMAV_SYSTEM_SOLVER == NUMAV_LDLT_SOLVER
            _a_diag.fill(_cmplx_t(0.0, 0.0));
        #endif

        // frequency calculations
        const double freq = _freq_steps[fi];
        if (freq == 0.0) {
            _x.fill(_cmplx_t(0.0, 0.0));

            // write solution to nmvr file
            nmvr::write_data_chunk_body(
                _nvmr_file,
                _ni_count * sizeof(_cmplx_t),
                _x.data()
            );
            
            // progress bar tick
            bar.set_progress(static_cast<size_t>(fi + 1UL));

            continue;
        }
        const double omega = 2.0 * std::numbers::pi * freq;
        const double omega_squared = omega * omega;

        // add point volume velocity to b vector
        for (size_t pvni = 0UL; pvni != _pvni_count; ++pvni)
        {
            const _cmplx_t volvel = (_pvni_to_forc_fi_part[pvni])(freq);
            *_pvni_to_ptr_in_b[pvni] += _cmplx_t(0.0, -omega) * volvel;
        }

        // add surface velocity to b vector
        for (size_t ispgv = 0UL; ispgv != _ispgv_count; ++ispgv)
        {
            const _cmplx_t velocity = (_ispgv_to_velocity[ispgv])(freq);
            const _cmplx_t fd_part = _cmplx_t(0.0, -omega) * velocity;
            for (size_t fipi = 0UL;
                fipi != _ispgv_to_ptr_in_b[ispgv].size(); ++fipi)
            {
                *_ispgv_to_ptr_in_b[ispgv][fipi] +=
                    _ispgv_to_forc_fi_part[ispgv][fipi] * fd_part;
            }
        }

        // add damping matrix to a
        for (size_t ispgi = 0UL; ispgi != _ispgi_count; ++ispgi) {
            const _cmplx_t impedance_value = _ispgi_to_impedance[ispgi](freq);
            const _cmplx_t damp_fd_part = _cmplx_t(0.0, omega)/impedance_value;
            
            for (size_t fipi = 0UL;
                fipi != _ispgi_to_ptr_in_a[ispgi].size(); ++fipi
            ) {
                *_ispgi_to_ptr_in_a[ispgi][fipi] +=
                    _ispgi_to_damp_fi_part[ispgi][fipi] * damp_fd_part;
            }
        }

        // add stiffness and mass matrix to a
        for (size_t ivpg = 0UL; ivpg != _ivpg_count; ++ivpg) {
            const _cmplx_t density_value =
                (_ivpg_to_volprop[ivpg].density)(freq);
            const _cmplx_t soundspeed_value =
                (_ivpg_to_volprop[ivpg].soundspeed)(freq);

            const _cmplx_t stif_fd_part = _cmplx_t(1.0, 0.0) / density_value;
            const _cmplx_t mass_fd_part = - omega_squared /
                (density_value * std::pow(soundspeed_value, 2UL));
            
            for (size_t fipi = 0UL;
                fipi != _ivpg_to_ptr_in_a[ivpg].size(); ++fipi) {
                *_ivpg_to_ptr_in_a[ivpg][fipi] +=
                    _ivpg_to_stif_fi_part[ivpg][fipi] * stif_fd_part +
                    _ivpg_to_mass_fi_part[ivpg][fipi] * mass_fd_part;
            }
        }

        // add pressure to a and b
        for (size_t pvi = 0UL; pvi != _pvi_count; ++pvi) {
            const _cmplx_t pressure = (_pvi_to_pressure[pvi])(freq);
            for (size_t fipi = 0UL;
                fipi != _pvi_to_ptr_in_a[pvi].size(); ++fipi
            ) {
                *_pvi_to_ptr_in_a[pvi][fipi] += PENALTY_METHOD_CONSTANT;
                *_pvi_to_ptr_in_b[pvi][fipi] +=
                    PENALTY_METHOD_CONSTANT * pressure;
            }
        }

        // solve
        #if NUMAV_SYSTEM_SOLVER == NUMAV_LDLT_SOLVER
            solve_using_ldlt_solver(
                _solver,
                _b_vals,
                _b_row_idx,
                _b_dense
            );
        #elif NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
            solve_using_eigen(
                _a_vals,
                _nnz_rowcol_idx_pairs,
                _b_vals,
                _b_row_idx,
                _ni_count,
                _x.data()
            );
        #elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
            solve_using_onemkl(
                _dss_handle,
                _a_vals,
                _b_vals,
                _b_row_idx,
                _b_dense,
                _x.data()
            );
        #else
            static_assert(false, "Invalid NUMAV_SYSTEM_SOLVER.");
        #endif

        // write solution to nmvr file
        nmvr::write_data_chunk_body(
            _nvmr_file,
            _ni_count * sizeof(_cmplx_t),
            _x.data()
        );
        
        // progress bar tick
        bar.set_progress(static_cast<size_t>(fi + 1UL));
    }
    #if NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
        constexpr MKL_INT options = NUMAV_MKL_OPTIONS;
        _INTEGER_t error_id = dss_delete(_dss_handle, options);
        if (error_id != MKL_DSS_SUCCESS) { print_dss_error(error_id); }
    #endif
    _end_nmvr_file();
    
    // end progress bar
    indicators::show_console_cursor(true);
    
    // print finish
    auto end_time = std::chrono::system_clock::now();
    auto time_t_end = std::chrono::system_clock::to_time_t(end_time);
    std::cout << "Solver ended at: " <<
        std::put_time(std::localtime(&time_t_end), "%Hh:%Mm:%Ss") << "\n";

    _did_run = true;
}

// explicit instantiation declarations
INSTANTIATE_SIMULATION_CLASS

} // namespace numav
