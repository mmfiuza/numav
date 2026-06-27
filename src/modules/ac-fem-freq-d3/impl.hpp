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
#include "indicators/progress_bar.hpp"
#include "H5Cpp.h"
#include "boost/bimap.hpp"

namespace numav {

template<ElementOrder O>
class SimulationAcFemFreqD3Tet<O>::Impl
{
public:
    Impl();
    ~Impl();
    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;
    Impl(Impl&&) noexcept;
    Impl& operator=(Impl&&) noexcept;

    NUMAV_SIM_AC_FEM_FREQ_D3_PUBLIC_METHODS

private:
    void _load_bdf(const char* const path);
    void _generate_extra_nodes();
    uint64_t _get_closest_point(const std::array<Float,3UL> point_coords);
    void _check_if_mesh_is_defined();
    void _check_if_did_run();
    void _validate_espg(const uint64_t espg);
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
    H5::DataSet _begin_hdf5_file();
    void _write_simulation_inputs_to_hdf5_file();
    void _write_solution_for_one_freq(H5::DataSet& ds, const uint64_t fi);
    void _clear_data_not_used_in_freq_iterations();
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

    enum class _FreqTypeDefinedByUser : uint64_t {
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

    H5::H5File _hdf5_file;

    std::unordered_set<uint64_t> _existing_evpg;
    std::unordered_set<uint64_t> _existing_espg;

    boost::bimap<uint64_t, uint64_t> _evpg_ivpg_bimap;
    boost::bimap<uint64_t, uint64_t> _espg_ispgi_bimap;
    boost::bimap<uint64_t, uint64_t> _espg_ispgv_bimap;
    boost::bimap<uint64_t, uint64_t> _espg_ispgp_bimap;

    std::vector<_VolProp> _ivpg_to_volprop;
    std::vector<FuncFloatToCmplx> _ispgi_to_impedance;
    std::vector<FuncFloatToCmplx> _ispgv_to_velocity;
    std::vector<FuncFloatToCmplx> _ispgp_to_pressure;

    std::string _hdf5_file_path;

    std::vector<uint64_t> _vpi_to_ni;
    std::vector<FuncFloatToCmplx> _vpi_to_volvel;
    std::vector<uint64_t> _ppi_to_ni;
    std::vector<FuncFloatToCmplx> _ppi_to_pressure;

    std::vector<std::array<Float, DIM>> _receiver_points;

    // members allocated during mesh load
    fz::SafePtr<std::array<Float, DIM>> _ni_to_coords;
    fz::SafePtr<std::array<uint64_t, ENIS_COUNT<O>>> _sei_to_ni;
    fz::SafePtr<std::array<uint64_t, ENIV_COUNT<O>>> _vei_to_ni;
    fz::SafePtr<uint64_t> _sei_to_espg;
    fz::SafePtr<uint64_t> _vei_to_evpg;

    // members allocated during the simulation run
    fz::SafePtr<uint64_t> _isei_to_sei;
    fz::SafePtr<uint64_t> _vsei_to_sei;
    fz::SafePtr<uint64_t> _psei_to_sei;
    fz::SafePtr<uint64_t> _vei_to_ivpg;
    fz::SafePtr<uint64_t> _isei_to_ispgi;
    fz::SafePtr<uint64_t> _vsei_to_ispgv;
    fz::SafePtr<uint64_t> _psei_to_ispgp;
    fz::SafePtr<std::pair<uint64_t, uint64_t>> _ni_connections;
    
    // members used during frequency iterations
    fz::SafePtr<Float> _fi_to_freq;
    fz::SafePtr<fz::SafePtr<Float>> _ivpg_to_stif_fi_part;
    fz::SafePtr<fz::SafePtr<Float>> _ivpg_to_mass_fi_part;
    fz::SafePtr<fz::SafePtr<Cmplx*>> _ivpg_to_ptr_in_a;
    fz::SafePtr<fz::SafePtr<Float>> _ispgi_to_damp_fi_part;
    fz::SafePtr<fz::SafePtr<Cmplx*>> _ispgi_to_ptr_in_a;
    fz::SafePtr<Cmplx*> _vpi_to_ptr_in_b;
    fz::SafePtr<fz::SafePtr<Float>> _ispgv_to_forc_fi_part;
    fz::SafePtr<fz::SafePtr<Cmplx*>> _ispgv_to_ptr_in_b;
    fz::SafePtr<FuncFloatToCmplx> _apvi_to_pressure;
    fz::SafePtr<fz::SafePtr<Cmplx*>> _apvi_to_ptr_in_a;
    fz::SafePtr<fz::SafePtr<Cmplx*>> _apvi_to_ptr_in_b;
    fz::SafePtr<Cmplx> _a_vals;
    fz::SafePtr<uint64_t> _b_row_idx;
    fz::SafePtr<Cmplx> _b_vals;
    fz::SafePtr<Cmplx> _x;

    uint64_t _fi_count;
    uint64_t _ni_count;
    uint64_t _sei_count;
    uint64_t _vei_count;
    uint64_t _isei_count;
    uint64_t _vsei_count;
    uint64_t _psei_count;
    uint64_t _ivpg_count;
    uint64_t _ispgi_count;
    uint64_t _ispgp_count;
    uint64_t _ispgv_count;
    uint64_t _vpi_count;
    uint64_t _ppi_count;
    uint64_t _apvi_count;
    uint64_t _ri_count;

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
