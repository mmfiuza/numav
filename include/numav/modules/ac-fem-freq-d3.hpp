// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include <unordered_map>
#include <unordered_set>

#include "Eigen/Eigen"

#define SAFE_PTR_DEBUG
#include "SafePtr.hpp"

namespace numav {

template<ElementOrder O> constexpr size_t NODES_IN_2D_ELEM = [] {
    if constexpr (O == ElementOrder::O1) { return 3; }
    if constexpr (O == ElementOrder::O2) { return 6; }
    return 0;
}();

template<ElementOrder O> constexpr size_t EXTRA_NODES_IN_2D_ELEM = [] {
    return NODES_IN_2D_ELEM<O> - NODES_IN_2D_ELEM<ElementOrder::O1>;
}();

template<ElementOrder O> constexpr size_t NODES_IN_3D_ELEM = [] {
    if constexpr (O == ElementOrder::O1) { return 4;  }
    if constexpr (O == ElementOrder::O2) { return 10; }
    return 0;
}();

template<ElementOrder O> constexpr size_t EXTRA_NODES_IN_3D_ELEM = [] {
    return NODES_IN_3D_ELEM<O> - NODES_IN_3D_ELEM<ElementOrder::O1>;
}();

template<> class Result<
    Phenomenon::ACOUSTIC,
    NumericalMethod::FEM,
    Domain::FREQUENCY,
    Dimension::D3
> {
public:
    Result();
    Result(const size_t&, const size_t&); // todo: make private
    ~Result();
    // todo: rule of 5
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> _data;
private:
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
    void set_frequency_range(
        const double&,
        const double&
    );
    void load_mesh(
        const char* const
    );
    void add_volume_material(
        const size_t&,
        const std::function<std::complex<double>(const double&)>&,
        const std::function<std::complex<double>(const double&)>&
    );
    void add_sound_source(
        const TypeOfSource&,
        const std::array<double,3>&,
        const PhysicalQuantity&,
        const std::function<std::complex<double>(const double&)>&
    );
    void add_sound_source(
        const TypeOfSource&,
        const size_t&,
        const PhysicalQuantity&,
        const std::function<std::complex<double>(const double&)>&
    );
    void add_surface_specific_acoustic_impedance(
        const size_t&,
        const std::function<std::complex<double>(const double&)>&
    );
    Result<
        Phenomenon::ACOUSTIC,
        NumericalMethod::FEM,
        Domain::FREQUENCY,
        Dimension::D3
    > run();

private:
    #include "numav/typedefs.hpp"

    // volume element properties
    struct _VolProp {
        _FuncRealToCmplx density;
        _FuncRealToCmplx soundspeed;
    };

    size_t _node_count() const;
    size_t _sfc_elem_count() const;
    size_t _vol_elem_count() const;
    size_t _ivpg_count() const;
    void _load_bdf(const char* const);
    void _generate_extra_nodes();
    _idx_t _get_closest_point(const std::array<double,3>&);
    void _check_if_mesh_is_defined();
    void _check_if_it_can_run();
    void _define_freq_vector();
    void _organize_physical_group_data();
    void _analyze_sparsity();
    void _assemble_freq_independent_parts();
    void _solve();

    bool _is_mesh_defined;
    bool _is_any_source_defined;
    bool _is_freq_range_defined;

    double _freq_min;
    double _freq_max;
    fz::SafePtr<double> _freq_steps;

    fz::SafePtr<std::array<double,3>> _node_coords;
    fz::SafePtr<std::array<_idx_t,NODES_IN_2D_ELEM<O>>> _sfc_elem_node_idx;
    fz::SafePtr<std::array<_idx_t,NODES_IN_3D_ELEM<O>>> _vol_elem_node_idx;

    std::vector<std::tuple<_idx_t,_FuncRealToCmplx>> _point_volvel;
    std::vector<std::tuple<_idx_t,_FuncRealToCmplx>> _point_pressure;

    std::unordered_set<_epg_t> _existing_espg;
    std::unordered_set<_epg_t> _existing_evpg;

    fz::SafePtr<_epg_t> _elem_idx_to_espg;
    fz::SafePtr<_epg_t> _elem_idx_to_evpg;
    
    std::unordered_map<_epg_t,_FuncRealToCmplx> _espg_to_volvel;
    std::unordered_map<_epg_t,_FuncRealToCmplx> _espg_to_pressure;
    std::unordered_map<_epg_t,_FuncRealToCmplx> _espg_to_impedance;
    std::unordered_map<_epg_t,_VolProp>         _evpg_to_volprop;

    std::unordered_map<_epg_t,_ipg_t> _espg_to_ispg;
    std::unordered_map<_epg_t,_ipg_t> _evpg_to_ivpg;

    fz::SafePtr<_ipg_t> _elem_idx_to_ispg;
    fz::SafePtr<_ipg_t> _elem_idx_to_ivpg;
    
    fz::SafePtr<_FuncRealToCmplx> _ispg_to_volvel;
    fz::SafePtr<_FuncRealToCmplx> _ispg_to_pressure;
    fz::SafePtr<_FuncRealToCmplx> _ispg_to_impedance;
    fz::SafePtr<_VolProp>         _ivpg_to_volprop;

    fz::SafePtr<std::pair<_idx_t,_idx_t>> _nnz_rowcol_idx_pairs;

    fz::SafePtr<fz::SafePtr<double>> _ivpg_to_btb_detj_w_vals;
    fz::SafePtr<fz::SafePtr<double>> _ivpg_to_nnt_detj_w_vals;
    fz::SafePtr<fz::SafePtr<std::complex<double>*>> _ivpg_to_ptr_in_a;

    fz::SafePtr<std::complex<double>> _a_vals;

    Result<
        Phenomenon::ACOUSTIC,
        NumericalMethod::FEM,
        Domain::FREQUENCY,
        Dimension::D3
    > _result;
};

// explicit instantiation declarations
template class Simulation<
    Phenomenon::ACOUSTIC,
    NumericalMethod::FEM,
    Domain::FREQUENCY,
    Dimension::D3,
    ElementOrder::O1
>;
template class Simulation<
    Phenomenon::ACOUSTIC,
    NumericalMethod::FEM,
    Domain::FREQUENCY,
    Dimension::D3,
    ElementOrder::O2
>;

} // namespace numav
