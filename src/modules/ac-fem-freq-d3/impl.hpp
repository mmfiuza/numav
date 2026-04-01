// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav/numav.hpp"
#include "common/aliases.hpp"
#include "common/debug-macros.hpp"
#include "modules/ac-fem-freq-d3/macros.hpp"
#include "modules/ac-fem-freq-d3/constants.hpp"

#include <unordered_map>
#include <unordered_set>

#include "Eigen/Eigen"
#if NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
    #include "mkl_dss.h"
    #include "mkl_types.h"
#endif
#include "SafePtr.hpp"

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

    void set_frequency_range(
        const double&,
        const double&
    );
    void load_mesh(
        const char* const
    );
    void add_volume_material(
        const size_t&,
        const std::function<_cmplx_t(const double&)>&,
        const std::function<_cmplx_t(const double&)>&
    );
    void add_volume_material(
        const size_t&,
        const char* const,
        const std::function<_cmplx_t(const double&)>&
    );
    void add_volume_material(
        const size_t&,
        const std::function<_cmplx_t(const double&)>&,
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
        const std::function<_cmplx_t(const double&)>&
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
        const std::function<_cmplx_t(const double&)>&
    );
    void add_sound_source(
        const TypeOfSource&,
        const size_t&,
        const PhysicalQuantity&,
        const char* const
    );
    void add_surface_specific_acoustic_impedance(
        const size_t&,
        const std::function<_cmplx_t(const double&)>&
    );
    void add_surface_specific_acoustic_impedance(
        const size_t&,
        const char* const
    );
    void run();
    void export_result(const char* const);

    // volume element properties
    struct _VolProp {
        _FuncRealToCmplx density;
        _FuncRealToCmplx soundspeed;
    };

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
    
    void _load_bdf(const char* const);
    void _generate_extra_nodes();
    size_t _get_closest_point(const std::array<double,3>&);
    void _check_if_mesh_is_defined();
    void _validate_espg(const size_t&);
    void _check_if_it_can_run();
    void _define_freq_vector();
    void _organize_volume_physical_group_data();
    void _organize_pressure_physical_group_data();
    void _organize_impedance_physical_group_data();
    void _organize_velocity_physical_group_data();
    void _organize_physical_group_data();
    void _allocate_a();
    void _allocate_b();
    void _assemble_fi_part_for_point_velocity();
    void _assemble_fi_part_for_sfc_velocity();
    void _assemble_fi_part_for_sfc_impedance();
    void _assemble_fi_part_for_vol_elements();
    void _assemble_fi_part_for_pressure();
    void _assemble_freq_independent_parts();
    void _solve_systems();

private:
    bool _is_mesh_defined;
    bool _is_any_source_defined;
    bool _is_freq_range_defined;
    bool _did_run;
    
    double _freq_min;
    double _freq_max;
    fz::SafePtr<double> _freq_steps;

    fz::SafePtr<std::array<double,DIM>> _ni_to_coords;
    fz::SafePtr<std::array<size_t,NODES_IN_SFC_ELEM<O>>> _sei_to_ni;
    fz::SafePtr<std::array<size_t,NODES_IN_VOL_ELEM<O>>> _vei_to_ni;

    std::vector<std::tuple<size_t,_FuncRealToCmplx>> _point_volvel;
    std::vector<std::tuple<size_t,_FuncRealToCmplx>> _point_pressure;

    std::unordered_set<size_t> _existing_evpg;
    std::unordered_set<size_t> _existing_espg;

    fz::SafePtr<size_t> _sei_to_espg;
    fz::SafePtr<size_t> _vei_to_evpg;
    
    std::unordered_map<size_t,_VolProp> _evpg_to_volprop;
    std::unordered_map<size_t,_FuncRealToCmplx> _espg_to_impedance;
    std::unordered_map<size_t,_FuncRealToCmplx> _espg_to_velocity;
    std::unordered_map<size_t,_FuncRealToCmplx> _espg_to_pressure;

    std::unordered_map<size_t,size_t> _espg_to_ispg;
    std::unordered_map<size_t,size_t> _evpg_to_ivpg;

    fz::SafePtr<size_t> _isei_to_sei;
    fz::SafePtr<size_t> _vsei_to_sei;
    fz::SafePtr<size_t> _psei_to_sei;

    fz::SafePtr<size_t> _vei_to_ivpg;
    fz::SafePtr<size_t> _isei_to_ispgi;
    fz::SafePtr<size_t> _vsei_to_ispgv;
    fz::SafePtr<size_t> _psei_to_ispgp;

    fz::SafePtr<std::pair<size_t,size_t>> _nnz_rowcol_idx_pairs;
    fz::SafePtr<_cmplx_t> _a_vals;
    fz::SafePtr<size_t> _b_row_idx;
    fz::SafePtr<_cmplx_t> _b_vals;
    
    fz::SafePtr<_VolProp> _ivpg_to_volprop;
    fz::SafePtr<_FuncRealToCmplx> _ispgi_to_impedance;
    fz::SafePtr<_FuncRealToCmplx> _ispgv_to_velocity;
    fz::SafePtr<_FuncRealToCmplx> _ispgp_to_pressure;

    fz::SafePtr<fz::SafePtr<double>> _ivpg_to_stif_fi_part;
    fz::SafePtr<fz::SafePtr<double>> _ivpg_to_mass_fi_part;
    fz::SafePtr<fz::SafePtr<_cmplx_t*>> _ivpg_to_ptr_in_a;

    fz::SafePtr<fz::SafePtr<double>> _ispgi_to_damp_fi_part;
    fz::SafePtr<fz::SafePtr<_cmplx_t*>> _ispgi_to_ptr_in_a;

    fz::SafePtr<_FuncRealToCmplx> _pvni_to_forc_fi_part;
    fz::SafePtr<_cmplx_t*> _pvni_to_ptr_in_b;

    fz::SafePtr<fz::SafePtr<double>> _ispgv_to_forc_fi_part;
    fz::SafePtr<fz::SafePtr<_cmplx_t*>> _ispgv_to_ptr_in_b;
    
    fz::SafePtr<_FuncRealToCmplx> _pvi_to_pressure;
    fz::SafePtr<fz::SafePtr<_cmplx_t*>> _pvi_to_ptr_in_a;
    fz::SafePtr<fz::SafePtr<_cmplx_t*>> _pvi_to_ptr_in_b;

    Eigen::Matrix<_cmplx_t,Eigen::Dynamic,Eigen::Dynamic> _cmplx_pressure_amp;

    #if NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
        _MKL_DSS_HANDLE_t _dss_handle;
        fz::SafePtr<_cmplx_t> _b_dense;
    #endif
};

} // namespace numav
