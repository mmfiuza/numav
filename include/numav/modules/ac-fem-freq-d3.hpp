// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include <unordered_map>
#include <unordered_set>

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
    ~Result();
    // todo: rule of 5
private:
    double _data;
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
    void _load_bdf(const char* const);
    void _generate_extra_nodes();
    _idx_t _get_closest_point(const std::array<double,3>&);
    void _check_if_mesh_is_defined();
    void _check_if_it_can_run();
    void _define_freq_vector();
    void _organize_physical_group_data();
    void _analize_sparsity();
    void _solve();

    bool _is_mesh_defined;
    bool _is_any_source_defined;
    bool _is_freq_range_defined;

    double _freq_min;
    double _freq_max;
    fz::SafePtr<double> _freq_vec;

    fz::SafePtr<std::array<double,3>> _node_coords;
    fz::SafePtr<std::array<_idx_t,NODES_IN_2D_ELEM<O>>> _sfc_elem_node_idx;
    fz::SafePtr<std::array<_idx_t,NODES_IN_3D_ELEM<O>>> _vol_elem_node_idx;

    std::vector<std::tuple<_idx_t,_FuncRealToCmplx>> _point_volvel;
    std::vector<std::tuple<_idx_t,_FuncRealToCmplx>> _point_pressure;

    std::unordered_set<_epg_t> _existing_epg_sfc;
    std::unordered_set<_epg_t> _existing_epg_vol;

    fz::SafePtr<_epg_t> _epg_sfc_elem;
    fz::SafePtr<_epg_t> _epg_vol_elem;
    
    std::unordered_map<_epg_t,_FuncRealToCmplx> _epg_to_sfc_volvel;
    std::unordered_map<_epg_t,_FuncRealToCmplx> _epg_to_sfc_pressure;
    std::unordered_map<_epg_t,_FuncRealToCmplx> _epg_to_sfc_impedance;
    std::unordered_map<_epg_t,_VolProp>         _epg_to_volprop;

    std::unordered_map<_epg_t,_ipg_t> _epg_to_ipg_sfc;
    std::unordered_map<_epg_t,_ipg_t> _epg_to_ipg_vol;

    fz::SafePtr<_ipg_t> _ipg_sfc_elem;
    fz::SafePtr<_ipg_t> _ipg_vol_elem;
    
    fz::SafePtr<_FuncRealToCmplx> _ipg_to_sfc_volvel;
    fz::SafePtr<_FuncRealToCmplx> _ipg_to_sfc_pressure;
    fz::SafePtr<_FuncRealToCmplx> _ipg_to_sfc_impedance;
    fz::SafePtr<_VolProp>         _ipg_to_volprop;

    fz::SafePtr<_idx_t> _nonzero_row_idx;
    fz::SafePtr<_idx_t> _nonzero_col_idx;
    fz::SafePtr<_idx_t> _btb_detj_w_vals;
    fz::SafePtr<_idx_t> _nnt_detj_w_vals;
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
