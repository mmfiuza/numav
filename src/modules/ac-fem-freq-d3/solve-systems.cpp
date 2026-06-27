// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/impl.hpp"

#include <cmath>

#include "common/utils.hpp"
#include "common/log.hpp"

namespace numav {

template<ElementOrder O>
void SimulationAcFemFreqD3Tet<O>::Impl::_clear_data_not_used_in_freq_iterations()
{
    // TODO: call shrink_to_fit in all std::vectors here  

    _existing_evpg.clear();
    _existing_espg.clear();

    _evpg_ivpg_bimap.clear();
    _espg_ispgi_bimap.clear();
    _espg_ispgv_bimap.clear();
    _espg_ispgp_bimap.clear();
    
    _receiver_points.clear();
    _ni_to_coords.free();
    _sei_to_ni.free();
    _vei_to_ni.free();
    _sei_to_espg.free();
    _vei_to_evpg.free();
    _isei_to_sei.free();
    _vsei_to_sei.free();
    _psei_to_sei.free();
    _vei_to_ivpg.free();
    _isei_to_ispgi.free();
    _vsei_to_ispgv.free();
    _psei_to_ispgp.free();
    _ni_connections.free();
}

template<ElementOrder O>
void SimulationAcFemFreqD3Tet<O>::Impl::_solve_systems()
{
    H5::DataSet pressure_data_set = _begin_hdf5_file();
    _write_simulation_inputs_to_hdf5_file();

    _clear_data_not_used_in_freq_iterations();

    log::print_start_time();

    // start progress bar
    const uint64_t bar_progress_max =
        (_fi_to_freq[0UL] == 0_F) ? _fi_count - 1UL : _fi_count;
    uint64_t bar_progress = 
        log::start_progress_bar(_progress_bar, bar_progress_max);

    for (uint64_t fi = 0UL; fi != _fi_count; ++fi)
    {
        // frequency calculations
        const Float freq = _fi_to_freq[fi];
        const Float omega = 2_F * PI * freq;
        const Float omega_squared = omega * omega;
        if (freq == 0_F) {
            _x.fill(Cmplx(0_F, 0_F));
            goto write_results_to_file_and_continue;
        }

        // reset values of A matrix and b vector
        _a_vals.fill(Cmplx(0_F, 0_F));
        #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
            _a_diag.fill(Cmplx(0_F, 0_F));
        #endif
        _b_vals.fill(Cmplx(0_F, 0_F));

        // add point volume velocity to b vector
        for (uint64_t vpi = 0UL; vpi != _vpi_count; ++vpi) {
            const Cmplx volvel = (_vpi_to_volvel[vpi])(freq);
            *_vpi_to_ptr_in_b[vpi] += Cmplx(0_F, -omega) * volvel;
        }

        // add surface velocity to b vector
        for (uint64_t ispgv = 0UL; ispgv != _ispgv_count; ++ispgv)
        {
            const Cmplx velocity = (_ispgv_to_velocity[ispgv])(freq);
            const Cmplx fd_part = Cmplx(0_F, -omega) * velocity;
            const uint64_t fipi_count = _ispgv_to_ptr_in_b[ispgv].size();
            for (uint64_t fipi = 0UL; fipi != fipi_count; ++fipi) {
                *_ispgv_to_ptr_in_b[ispgv][fipi] +=
                    _ispgv_to_forc_fi_part[ispgv][fipi] * fd_part;
            }
        }

        // add damping matrix to a
        for (uint64_t ispgi = 0UL; ispgi != _ispgi_count; ++ispgi)
        {
            const Cmplx impedance = _ispgi_to_impedance[ispgi](freq);
            const Cmplx damp_fd_part = Cmplx(0_F, omega) / impedance;
            const uint64_t fipi_count = _ispgi_to_ptr_in_a[ispgi].size();
            for (uint64_t fipi = 0UL; fipi != fipi_count; ++fipi) {
                *_ispgi_to_ptr_in_a[ispgi][fipi] +=
                    _ispgi_to_damp_fi_part[ispgi][fipi] * damp_fd_part;
            }
        }

        // add stiffness and mass matrix to a
        for (uint64_t ivpg = 0UL; ivpg != _ivpg_count; ++ivpg)
        {
            const Cmplx density = (_ivpg_to_volprop[ivpg].density)(freq);
            const Cmplx soundspeed = (_ivpg_to_volprop[ivpg].soundspeed)(freq);
            const Cmplx stif_fd_part = Cmplx(1_F, 0_F) / density;
            const Cmplx mass_fd_part =
                - omega_squared / (density * soundspeed * soundspeed);
            
            const uint64_t fipi_count = _ivpg_to_ptr_in_a[ivpg].size();
            for (uint64_t fipi = 0UL; fipi != fipi_count; ++fipi) {
                *_ivpg_to_ptr_in_a[ivpg][fipi] +=
                    _ivpg_to_stif_fi_part[ivpg][fipi] * stif_fd_part +
                    _ivpg_to_mass_fi_part[ivpg][fipi] * mass_fd_part;
            }
        }

        // add pressure to a and b
        for (uint64_t apvi = 0UL; apvi != _apvi_count; ++apvi)
        {
            const Cmplx pressure = (_apvi_to_pressure[apvi])(freq);
            const uint64_t fipi_count = _apvi_to_ptr_in_a[apvi].size();
            for (uint64_t fipi = 0UL; fipi != fipi_count; ++fipi) {
                *_apvi_to_ptr_in_a[apvi][fipi] += PENALTY_METHOD_CONSTANT;
                *_apvi_to_ptr_in_b[apvi][fipi] +=
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

        write_results_to_file_and_continue:
        
        // write solution to hdf5 file
        _write_solution_for_one_freq(pressure_data_set, fi);
    }

    // terminate solver
    #if NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
        _terminate_onemkl_solver();
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_MUMPS
        _terminate_mumps_solver();
    #endif

    _hdf5_file.close();
    log::finish_progress_bar();
    log::print_finish_time();

    _did_run = true;
}

} // namespace numav

NUMAV_INSTANTIATE_SIM_AC_FEM_FREQ_D3
