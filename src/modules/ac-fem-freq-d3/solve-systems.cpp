// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/impl.hpp"

#include <cmath>
#include <numbers>

#include "common/nmvr-format.hpp"
#include "common/utils.hpp"
#include "common/log.hpp"

namespace numav {

template<ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_solve_systems()
{
    _begin_nmvr_file();
    _write_simulation_inputs_to_nmvr_file();
    nmvr::write_data_chunk_header(
        _nvmr_file,
        nmvr::COMPLEX_AMPLITUDE_OF_SOUND_PRESSURE_SOLUTION_CHUNK_ID,
        _ni_count * _freq_count * sizeof(std::complex<double>)
    );

    // print start time
    log::print_start_time();

    // start progress bar
    const size_t bar_progress_max =
        (_freq_steps[0UL] == 0_F) ? _freq_count - 1UL : _freq_count;
    size_t bar_progress = 
        log::start_progress_bar(_progress_bar, bar_progress_max);

    for (size_t fi = 0UL; fi != _freq_count; ++fi)
    {
        // frequency calculations
        const Float freq = _freq_steps[fi];
        const Float omega = 2_F * std::numbers::pi * freq;
        const Float omega_squared = omega * omega;
        if (freq == 0_F) {
            _x.fill(Cmplx(0_F, 0_F));
            goto write_results_to_nmvr_and_continue;
        }

        // reset values of A matrix and b vector
        _a_vals.fill(Cmplx(0_F, 0_F));
        #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
            _a_diag.fill(Cmplx(0_F, 0_F));
        #endif
        _b_vals.fill(Cmplx(0_F, 0_F));

        // add point volume velocity to b vector
        for (size_t pvni = 0UL; pvni != _pvni_count; ++pvni)
        {
            const Cmplx volvel = (_pvni_to_forc_fi_part[pvni])(freq);
            *_pvni_to_ptr_in_b[pvni] += Cmplx(0_F, -omega) * volvel;
        }

        // add surface velocity to b vector
        for (size_t ispgv = 0UL; ispgv != _ispgv_count; ++ispgv)
        {
            const Cmplx velocity = (_ispgv_to_velocity[ispgv])(freq);
            const Cmplx fd_part = Cmplx(0_F, -omega) * velocity;
            for (size_t fipi = 0UL;
                fipi != _ispgv_to_ptr_in_b[ispgv].size(); ++fipi)
            {
                *_ispgv_to_ptr_in_b[ispgv][fipi] +=
                    _ispgv_to_forc_fi_part[ispgv][fipi] * fd_part;
            }
        }

        // add damping matrix to a
        for (size_t ispgi = 0UL; ispgi != _ispgi_count; ++ispgi)
        {
            const Cmplx impedance_value = _ispgi_to_impedance[ispgi](freq);
            const Cmplx damp_fd_part = Cmplx(0_F, omega) / impedance_value;
            for (size_t fipi = 0UL;
                fipi != _ispgi_to_ptr_in_a[ispgi].size(); ++fipi
            ) {
                *_ispgi_to_ptr_in_a[ispgi][fipi] +=
                    _ispgi_to_damp_fi_part[ispgi][fipi] * damp_fd_part;
            }
        }

        // add stiffness and mass matrix to a
        for (size_t ivpg = 0UL; ivpg != _ivpg_count; ++ivpg)
        {
            const Cmplx density_value =
                (_ivpg_to_volprop[ivpg].density)(freq);
            const Cmplx soundspeed_value =
                (_ivpg_to_volprop[ivpg].soundspeed)(freq);

            const Cmplx stif_fd_part = Cmplx(1_F, 0_F) / density_value;
            const Cmplx mass_fd_part = - omega_squared / 
                (density_value * soundspeed_value * soundspeed_value);
            
            for (size_t fipi = 0UL;
                fipi != _ivpg_to_ptr_in_a[ivpg].size(); ++fipi
            ) {
                *_ivpg_to_ptr_in_a[ivpg][fipi] +=
                    _ivpg_to_stif_fi_part[ivpg][fipi] * stif_fd_part +
                    _ivpg_to_mass_fi_part[ivpg][fipi] * mass_fd_part;
            }
        }

        // add pressure to a and b
        for (size_t pvi = 0UL; pvi != _pvi_count; ++pvi)
        {
            const Cmplx pressure = (_pvi_to_pressure[pvi])(freq);
            for (size_t fipi = 0UL;
                fipi != _pvi_to_ptr_in_a[pvi].size(); ++fipi
            ) {
                *_pvi_to_ptr_in_a[pvi][fipi] += PENALTY_METHOD_CONSTANT;
                *_pvi_to_ptr_in_b[pvi][fipi] +=
                    PENALTY_METHOD_CONSTANT * pressure;
            }
        }

        // solve
        #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
            _solve_using_internal_solver();
        #elif NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
            _solve_using_eigen_solver();
        #elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
            _solve_using_onemkl_solver();
        #elif NUMAV_SYSTEM_SOLVER == NUMAV_MUMPS
            _solve_using_mumps_solver();
        #else
            static_assert(false, "Invalid NUMAV_SYSTEM_SOLVER.");
        #endif

        // increment progress bar
        log::increment_progress_bar(_progress_bar, bar_progress);

        write_results_to_nmvr_and_continue:

        // write solution to nmvr file
        nmvr::write_data_chunk_body(
            _nvmr_file,
            _ni_count * sizeof(std::complex<double>),
            static_cast_contiguous_data<std::complex<double>>(
                _x.data(), _ni_count
            ).get()
        );
    }
    _end_nmvr_file();

    // terminate solver
    #if NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
        _terminate_onemkl_solver();
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_MUMPS
        _terminate_mumps_solver();
    #endif
    
    // finish progress bar
    log::finish_progress_bar();

    // print finish
    log::print_finish_time();

    _did_run = true;
}

// explicit instantiation declarations
INSTANTIATE_SIMULATION_CLASS

} // namespace numav
