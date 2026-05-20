// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/numav.hpp"
#include "numav/aliases.hpp"
#include "modules/ac-fem-freq-d3/macros.hpp"
#include "modules/ac-fem-freq-d3/compile-options.hpp"
#include "modules/ac-fem-freq-d3/constants.hpp"

#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <memory>
#include <optional>

#include "Eigen/Eigen"
#include "Eigen/OrderingMethods"
#include "SafePtr.hpp"
#include <indicators/progress_bar.hpp>

namespace numav {

template<ElementOrder O>
class SimulationAcFemFreqD3<O>::Impl
{
public:
    Impl();
    ~Impl();
    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;
    Impl(Impl&&) noexcept;
    Impl& operator=(Impl&&) noexcept;

    void set_maximum_frequency(
        const Float freq_max
    );
    void set_frequency_range(
        const Float freq_min,
        const Float freq_max
    );
    void set_frequency_steps_count(
        const size_t freq_steps_count
    );
    void set_frequency_sampling_density(
        const FrequencySamplingDensity freq_sampling_density
    );
    void set_frequency_steps(
        const std::vector<Float>& freq_steps
    );
    void load_mesh(
        const char* const path_to_mesh
    );
    void add_volume_material(
        const size_t evpg,
        const FuncFloatToCmplx& density_func,
        const FuncFloatToCmplx& soundspeed_func
    );
    void add_volume_material(
        const size_t evpg,
        const char* const density_table,
        const FuncFloatToCmplx& soundspeed_func
    );
    void add_volume_material(
        const size_t evpg,
        const FuncFloatToCmplx& density_func,
        const char* const soundspeed_table
    );
    void add_volume_material(
        const size_t evpg,
        const char* const density_table,
        const char* const soundspeed_table
    );
    void add_sound_source(
        const TypeOfSource source_type,
        const std::array<Float,3UL> point_coords,
        const PhysicalQuantity pq_type,
        const FuncFloatToCmplx& pq_func
    );
    void add_sound_source(
        const TypeOfSource source_type,
        const std::array<Float,3UL> point_coords,
        const PhysicalQuantity pq_type,
        const char* const pq_table
    );
    void add_sound_source(
        const TypeOfSource source_type,
        const size_t espg,
        const PhysicalQuantity pq_type,
        const FuncFloatToCmplx& pq_func
    );
    void add_sound_source(
        const TypeOfSource source_type,
        const size_t espg,
        const PhysicalQuantity pq_type,
        const char* const pq_table
    );
    void add_receiver(
        const std::array<Float,DIM> point_coords
    );
    void add_surface_material(
        const size_t espg,
        const PhysicalQuantity pq_type,
        const FuncFloatToCmplx& pq_func
    );
    void add_surface_material(
        const size_t espg,
        const PhysicalQuantity pq_type,
        const char* const pq_table
    );
    void set_result_export_path(const char* const path_to_result);
    void run();

private:

    void _load_bdf(const char* const path);
    void _generate_extra_nodes();
    size_t _get_closest_point(const std::array<Float,3UL> point_coords);
    void _check_if_mesh_is_defined();
    void _validate_espg(const size_t espg);
    void _check_if_it_can_run();
    void _define_freq_vector();
    void _organize_volume_physical_group_data();
    void _organize_pressure_physical_group_data();
    void _organize_impedance_physical_group_data();
    void _organize_velocity_physical_group_data();
    void _organize_physical_group_data();
    void _allocate_a();
    void _allocate_b();
    void _allocate_x();
    void _assemble_fi_part_for_point_velocity();
    void _assemble_fi_part_for_sfc_velocity();
    void _assemble_fi_part_for_sfc_impedance();
    void _assemble_fi_part_for_vol_elements();
    void _assemble_fi_part_for_pressure();
    void _assemble_freq_independent_parts();
    void _solve_systems();
    void _post_process();
    void _begin_nmvr_file();
    void _write_simulation_inputs_to_nmvr_file();
    void _end_nmvr_file();
    #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
        void _define_sparsity_pattern_using_internal_solver();
        void _solve_using_internal_solver();
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
        void _define_sparsity_pattern_using_eigen_solver();
        void _solve_using_eigen_solver();
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
        void _define_sparsity_pattern_using_onemkl_solver();
        void _solve_using_onemkl_solver();
        void _terminate_onemkl_solver();
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_MUMPS
        void _define_sparsity_pattern_using_mumps_solver();
        void _solve_using_mumps_solver();
        void _terminate_mumps_solver();
    #endif

    enum class _FreqTypeDefinedByUser : size_t {
        UNDEFINED,
        MAXIMUM,
        RANGE,
        STEPS
    };
    
    // volume element properties
    struct _VolProp {
        FuncFloatToCmplx density;
        FuncFloatToCmplx soundspeed;
    };

    std::ofstream _nvmr_file;

    std::unordered_set<size_t> _existing_evpg;
    std::unordered_set<size_t> _existing_espg;
    std::unordered_map<size_t, _VolProp> _evpg_to_volprop;
    std::unordered_map<size_t, FuncFloatToCmplx> _espg_to_impedance;
    std::unordered_map<size_t, FuncFloatToCmplx> _espg_to_velocity;
    std::unordered_map<size_t, FuncFloatToCmplx> _espg_to_pressure;
    std::unordered_map<size_t, size_t> _espg_to_ispg;
    std::unordered_map<size_t, size_t> _evpg_to_ivpg;

    std::string _nmvr_file_path;

    std::vector<std::tuple<size_t, FuncFloatToCmplx>> _point_volvel;
    std::vector<std::tuple<size_t, FuncFloatToCmplx>> _point_pressure;
    std::vector<std::array<Float, DIM>> _receiver_points;

    fz::SafePtr<Float> _freq_steps;
    fz::SafePtr<std::array<Float, DIM>> _ni_to_coords;
    fz::SafePtr<std::array<size_t, NODES_IN_SFC_ELEM<O>>> _sei_to_ni;
    fz::SafePtr<std::array<size_t, NODES_IN_VOL_ELEM<O>>> _vei_to_ni;
    fz::SafePtr<size_t> _sei_to_espg;
    fz::SafePtr<size_t> _vei_to_evpg;
    fz::SafePtr<size_t> _isei_to_sei;
    fz::SafePtr<size_t> _vsei_to_sei;
    fz::SafePtr<size_t> _psei_to_sei;
    fz::SafePtr<size_t> _vei_to_ivpg;
    fz::SafePtr<size_t> _isei_to_ispgi;
    fz::SafePtr<size_t> _vsei_to_ispgv;
    fz::SafePtr<size_t> _psei_to_ispgp;
    fz::SafePtr<_VolProp> _ivpg_to_volprop;
    fz::SafePtr<FuncFloatToCmplx> _ispgi_to_impedance;
    fz::SafePtr<FuncFloatToCmplx> _ispgv_to_velocity;
    fz::SafePtr<FuncFloatToCmplx> _ispgp_to_pressure;
    fz::SafePtr<fz::SafePtr<Float>> _ivpg_to_stif_fi_part;
    fz::SafePtr<fz::SafePtr<Float>> _ivpg_to_mass_fi_part;
    fz::SafePtr<fz::SafePtr<Cmplx*>> _ivpg_to_ptr_in_a;
    fz::SafePtr<fz::SafePtr<Float>> _ispgi_to_damp_fi_part;
    fz::SafePtr<fz::SafePtr<Cmplx*>> _ispgi_to_ptr_in_a;
    fz::SafePtr<FuncFloatToCmplx> _pvni_to_forc_fi_part;
    fz::SafePtr<Cmplx*> _pvni_to_ptr_in_b;
    fz::SafePtr<fz::SafePtr<Float>> _ispgv_to_forc_fi_part;
    fz::SafePtr<fz::SafePtr<Cmplx*>> _ispgv_to_ptr_in_b;
    fz::SafePtr<FuncFloatToCmplx> _pvi_to_pressure;
    fz::SafePtr<fz::SafePtr<Cmplx*>> _pvi_to_ptr_in_a;
    fz::SafePtr<fz::SafePtr<Cmplx*>> _pvi_to_ptr_in_b;
    fz::SafePtr<std::pair<size_t, size_t>> _ni_connections;
    fz::SafePtr<Cmplx> _a_vals;
    fz::SafePtr<size_t> _b_row_idx;
    fz::SafePtr<Cmplx> _b_vals;
    fz::SafePtr<Cmplx> _x;

    size_t _ni_count;
    size_t _sei_count;
    size_t _vei_count;
    size_t _isei_count;
    size_t _vsei_count;
    size_t _psei_count;
    size_t _ivpg_count;
    size_t _ispgi_count;
    size_t _pvni_count;
    size_t _ispgv_count;
    size_t _ppni_count;
    size_t _ispgp_count;
    size_t _pvi_count;
    size_t _freq_count;
    size_t _ri_count;

    Float _freq_min;
    Float _freq_max;
    FrequencySamplingDensity _freq_sampling_density;
    _FreqTypeDefinedByUser _freq_type_defined_by_user;

    std::unique_ptr<indicators::ProgressBar> _progress_bar;

    bool _is_mesh_defined;
    bool _is_any_source_defined;
    bool _did_run;

    #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
        LdltSolver<Cmplx> _solver;
        fz::SafePtr<Cmplx> _b_dense;
        fz::SafePtr<Cmplx> _a_diag;
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
        std::optional<Eigen::Map<
            Eigen::SparseMatrix<Cmplx, Eigen::ColMajor, Eigen::Index>
        >> _a_eigen;
        std::array<Eigen::Index, 2UL> _b_col_idx_signed;
        std::optional<Eigen::Map<
            Eigen::SparseMatrix<Cmplx, Eigen::ColMajor, Eigen::Index>
        >> _b_eigen;
        std::optional<Eigen::Map<
            Eigen::Matrix<Cmplx, Eigen::Dynamic, 1UL>
        >> _x_eigen;
        fz::SafePtr<Eigen::Index> _a_row_idx;
        fz::SafePtr<Eigen::Index> _a_col_idx;
        fz::SafePtr<Eigen::Index> _b_row_idx_signed;
        std::unique_ptr<
            Eigen::SparseLU<
                Eigen::SparseMatrix<Cmplx, Eigen::ColMajor, Eigen::Index>,
                Eigen::COLAMDOrdering<Eigen::Index>
            >
        > _solver;
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
        fz::SafePtr<Cmplx> _b_dense;
        _MKL_DSS_HANDLE_t _dss_handle;
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_MUMPS
        ZMUMPS_STRUC_C _solver;
        fz::SafePtr<MUMPS_INT> _a_row_idx;
        fz::SafePtr<MUMPS_INT> _a_col_idx;
        fz::SafePtr<Cmplx> _b_dense;
    #endif
};

} // namespace numav
