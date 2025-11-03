// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/impl.hpp"

#include <tuple>
#include <cmath>
#include <numbers>
#include <algorithm>
#include <limits>
#include <chrono>
#include <iomanip>
#include <set>

#include "common/exception.hpp"
#include "common/log.hpp"
#include "common/hash-functions.hpp"
#include "common/maths.hpp"
#include "common/utils.hpp"

namespace numav {

// constants
static constexpr size_t DIM = DIM_COUNT<Dimension::D3>;
static constexpr double AREA_REF_TRIG = 1.0 / 2.0;
static constexpr double VOLUME_REF_TET = 1.0 / 6.0;
static constexpr double PENALTY_METHOD_CONSTANT = 1e12;

template<ElementOrder O> constexpr size_t NGP_FORC = [] {
    if constexpr (O == ElementOrder::O1) { return 1; }
    if constexpr (O == ElementOrder::O2) { return 3; }
}();

template<ElementOrder O> constexpr size_t NGP_DAMP = [] {
    if constexpr (O == ElementOrder::O1) { return 1; }
    if constexpr (O == ElementOrder::O2) { return 3; }
}();

template<ElementOrder O> constexpr size_t NGP_STIF = [] {
    if constexpr (O == ElementOrder::O1) { return 1; }
    if constexpr (O == ElementOrder::O2) { return 4; }
}();

template<ElementOrder O> constexpr size_t NGP_MASS = [] {
    if constexpr (O == ElementOrder::O1) { return 4;  }
    if constexpr (O == ElementOrder::O2) { return 15; }
}();

template<size_t N>
constexpr std::array<std::array<double,2>,N> GAUSS_POINTS_SFC = [] {
    if constexpr (N == 1) {
        constexpr double a = 1.0 / 3.0;
        return std::array<std::array<double,2>,N>({ {a, a} });
    }
    if constexpr (N == 3) {
        constexpr double a = 1.0 / 6.0;
        constexpr double b = 2.0 / 3.0;
        return std::array<std::array<double,2>,N>{{
            {a,a}, {b,a}, {a,b}
        }};
    }
}();

template<size_t N>
constexpr std::array<std::array<double,DIM>,N> GAUSS_POINTS_VOL = [] {
    if constexpr (N == 1) {
        constexpr double a = 1.0 / 4.0;
        return std::array<std::array<double,DIM>,N>({ {a, a, a} });
    }
    if constexpr (N == 4) {
        constexpr double a = (5.0 -     std::sqrt(5.0)) / 20.0;
        constexpr double b = (5.0 + 3.0*std::sqrt(5.0)) / 20.0;
        return std::array<std::array<double,DIM>,N>{{
            {a,a,a}, {b,a,a}, {a,b,a}, {a,a,b}
        }};
    }
    if constexpr (N == 5) {
        constexpr double a = 1.0 / 4.0;
        constexpr double b = 1.0 / 6.0;
        constexpr double c = 1.0 / 2.0;
        return std::array<std::array<double,DIM>,N>{{
            {a,a,a},
            {b,b,b}, {c,b,b}, {b,c,b}, {b,b,c}
        }};
    }
    if constexpr (N == 11) {
        constexpr double a =  1.0 / 4.0;
        constexpr double b =  1.0 / 14.0;
        constexpr double c = 11.0 / 14.0;
        constexpr double d = 0.399403576166799;
        constexpr double e = 0.100596423833201;
        return std::array<std::array<double,DIM>,N>{{
            {a,a,a},
            {b,b,b}, {c,b,b}, {b,c,b}, {b,b,c},
            {d,d,e}, {d,e,d}, {e,d,d}, {d,e,e}, {e,d,e}, {e,e,d}
        }};
    }
    if constexpr (N == 15) {
        constexpr double a = 1.0 / 4.0;
        constexpr double b = ( 7.0 +     std::sqrt(15.0)) / 34.0;
        constexpr double c = ( 7.0 -     std::sqrt(15.0)) / 34.0;
        constexpr double d = (13.0 - 3.0*std::sqrt(15.0)) / 34.0;
        constexpr double e = (13.0 + 3.0*std::sqrt(15.0)) / 34.0;
        constexpr double f = ( 5.0 -     std::sqrt(15.0)) / 20.0;
        constexpr double g = ( 5.0 +     std::sqrt(15.0)) / 20.0;
        return std::array<std::array<double,DIM>,N>{{
            {a,a,a},
            {b,b,b}, {b,b,d}, {b,d,b}, {d,b,b}, 
            {c,c,c}, {c,c,e}, {c,e,c}, {e,c,c},
            {f,f,g}, {f,g,f}, {g,f,f}, {f,g,g}, {g,f,g}, {g,g,f}
        }};
    }
}();

template<size_t N>
constexpr std::array<double,N> GAUSS_WEIGHTS_SFC = [] {
    if constexpr (N == 1) {
        constexpr double a = 1.0 / 2.0;
        return std::array<double,N>({a});
    }
    if constexpr (N == 3) {
        constexpr double a = 1.0 / 6.0;
        return std::array<double,N>({a,a,a});
    }
}();

template<size_t N>
constexpr std::array<double,N> GAUSS_WEIGHTS_VOL = [] {
    if constexpr (N == 1) {
        constexpr double a = 1.0 / 6.0;
        return std::array<double,N>({a});
    }
    if constexpr (N == 4) {
        constexpr double a = 1.0 / 24.0;
        return std::array<double,N>({a,a,a,a});
    }
    if constexpr (N == 5) {
        constexpr double a = -2.0 / 15.0;
        constexpr double b = 3.0 / 40.0;
        return std::array<double,N>({a,b,b,b,b});
    }
    if constexpr (N == 11) {
        constexpr double a = -74.0 / 5625.0;
        constexpr double b = 343.0 / 45000.0;
        constexpr double c = 28.0 / 1125.0;
        return std::array<double,N>({a,b,b,b,b,c,c,c,c,c,c});
    }
    if constexpr (N == 15) {
        constexpr double a = 8.0 / 405.0;
        constexpr double b = (2665.0 - 14.0*std::sqrt(15.0)) / 226800.0;
        constexpr double c = (2665.0 + 14.0*std::sqrt(15.0)) / 226800.0;
        constexpr double d = 5.0 / 567.0;
        return std::array<double,N>({a,b,b,b,b,c,c,c,c,d,d,d,d,d,d});
    }
}();

// general shape_func_sfc declaration for all orders
template<ElementOrder O>
Eigen::Matrix<double, NODES_IN_SFC_ELEM<O>, 1> shape_func_sfc(
    const double&, const double&
);

template<>
Eigen::Matrix<double, NODES_IN_SFC_ELEM<ElementOrder::O1>, 1>
shape_func_sfc<ElementOrder::O1>(
    const double& xi0, const double& xi1
) {
    const double xi2 = 1 - xi0 - xi1;
    return Eigen::Matrix<double, NODES_IN_SFC_ELEM<ElementOrder::O1>, 1>(
        xi0,
        xi1,
        xi2
    );
}

template<>
Eigen::Matrix<double, NODES_IN_SFC_ELEM<ElementOrder::O2>, 1>
shape_func_sfc<ElementOrder::O2>(
    const double& xi0, const double& xi1
) {
    const double xi2 = 1 - xi0 - xi1;
    return Eigen::Matrix<double, NODES_IN_SFC_ELEM<ElementOrder::O2>, 1>(
        xi0*(2*xi0-1),
        xi1*(2*xi1-1),
        xi2*(2*xi2-1),
        4*xi0*xi1,
        4*xi0*xi2,
        4*xi1*xi2
    );
}

// general shape_func_sfc_gradient declaration for all orders
template<ElementOrder O>
Eigen::Matrix<double, 2, NODES_IN_SFC_ELEM<O>> shape_func_sfc_gradient(
    const double&, const double&
);

template<>
Eigen::Matrix<double, 2, NODES_IN_SFC_ELEM<ElementOrder::O1>>
shape_func_sfc_gradient<ElementOrder::O1>(
    const double& xi0, const double& xi1
) {
    (void)xi0;
    (void)xi1;
    return Eigen::Matrix<double, 2, NODES_IN_SFC_ELEM<ElementOrder::O1>> {
        {1, 0, -1},
        {0, 1, -1}
    };
}

template<>
Eigen::Matrix<double, 2, NODES_IN_SFC_ELEM<ElementOrder::O2>>
shape_func_sfc_gradient<ElementOrder::O2>(
    const double& xi0, const double& xi1
) {
    const double xi2 = 1 - xi0 - xi1;
    return Eigen::Matrix<double, 2, NODES_IN_SFC_ELEM<ElementOrder::O2>> {
        {4*xi0-1, 0      , 1-4*xi2, 4*xi1, 4*(xi2-xi0), -4*xi1     },
        {0      , 4*xi1-1, 1-4*xi2, 4*xi0, -4*xi0     , 4*(xi2-xi1)}
    };
}

// general shape_func_vol declaration for all orders
template<ElementOrder O>
Eigen::Matrix<double, NODES_IN_VOL_ELEM<O>, 1> shape_func_vol(
    const double&, const double&, const double&
);

template<>
Eigen::Matrix<double, NODES_IN_VOL_ELEM<ElementOrder::O1>, 1>
shape_func_vol<ElementOrder::O1>(
    const double& xi0, const double& xi1, const double& xi2
) {
    const double xi3 = 1 - xi0 - xi1 - xi2;
    return Eigen::Matrix<double, NODES_IN_VOL_ELEM<ElementOrder::O1>, 1>(
        xi0,
        xi1,
        xi2,
        xi3
    );
}

template<>
Eigen::Matrix<double, NODES_IN_VOL_ELEM<ElementOrder::O2>, 1>
shape_func_vol<ElementOrder::O2>(
    const double& xi0, const double& xi1, const double& xi2
) {
    const double xi3 = 1 - xi0 - xi1 - xi2;
    return Eigen::Matrix<double, NODES_IN_VOL_ELEM<ElementOrder::O2>, 1>(
        xi0*(2*xi0-1),
        xi1*(2*xi1-1),
        xi2*(2*xi2-1),
        xi3*(2*xi3-1),
        4*xi0*xi1,
        4*xi0*xi2,
        4*xi0*xi3,
        4*xi1*xi2,
        4*xi1*xi3,
        4*xi2*xi3
    );
}

// general shape_func_vol_gradient declaration for all orders
template<ElementOrder O>
Eigen::Matrix<double, DIM, NODES_IN_VOL_ELEM<O>> shape_func_vol_gradient(
    const double&, const double&, const double&
);

template<>
Eigen::Matrix<double, DIM, NODES_IN_VOL_ELEM<ElementOrder::O1>>
shape_func_vol_gradient<ElementOrder::O1>(
    const double& xi0, const double& xi1, const double& xi2
) {
    (void)xi0;
    (void)xi1;
    (void)xi2;
    return Eigen::Matrix<double, DIM, NODES_IN_VOL_ELEM<ElementOrder::O1>> {
        {1, 0, 0, -1},
        {0, 1, 0, -1},
        {0, 0, 1, -1}
    };
}

template<>
Eigen::Matrix<double, DIM, NODES_IN_VOL_ELEM<ElementOrder::O2>>
shape_func_vol_gradient<ElementOrder::O2>(
    const double& xi0, const double& xi1, const double& xi2
) {
    const double xi3 = 1 - xi0 - xi1 - xi2;
    return Eigen::Matrix<double, DIM, NODES_IN_VOL_ELEM<ElementOrder::O2>> {
        {4*xi0-1, 0, 0, 1-4*xi3, 4*xi1, 4*xi2, 4*(xi3-xi0), 0, -4*xi1, -4*xi2},
        {0, 4*xi1-1, 0, 1-4*xi3, 4*xi0, 0, -4*xi0, 4*xi2, 4*(xi3-xi1), -4*xi2},
        {0, 0, 4*xi2-1, 1-4*xi3, 0, 4*xi0, -4*xi0, 4*xi1, -4*xi1, 4*(xi3-xi2)}
    };
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_check_if_it_can_run() {
    _check_if_mesh_is_defined();
    if (!_is_any_source_defined){
        error("No sound source was defined."
            " Call add_sound_source to do so.");
    }
    if (!_is_freq_range_defined){
        error("Simulation frequency range was not defined."
            " Call set_frequency_range to do so.");
    }
    for (auto& evpg : _existing_evpg) {
        if (!_evpg_to_volprop.contains(evpg)) {
            error("Volume tag {} was not assigned."
            " Call add_volume_material to do so.", evpg);
        }
    }
    if (_did_run) {
        error("This Simulation has already been run.");
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_define_freq_vector() {
    // todo: decide number here
    // todo: make it not linear
    _freq_count = 1000;
    _freq_steps = linspace(_freq_min, _freq_max, _freq_count);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_organize_volume_physical_group_data()
{
    _ivpg_to_volprop = fz::SafePtr<_VolProp>(_ivpg_count);
    for (const auto& [evpg, volprop] : _evpg_to_volprop) {
        const size_t ivpg = _evpg_to_ivpg.at(evpg);
        _ivpg_to_volprop[ivpg] = volprop;
    }
    _vei_to_ivpg = fz::SafePtr<size_t>(_vei_count);
    for (size_t vei=0; vei!=_vei_count; ++vei) {
        _vei_to_ivpg[vei] = _evpg_to_ivpg.at(_vei_to_evpg[vei]);
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_organize_impedance_physical_group_data()
{
    _ispgi_to_impedance = fz::SafePtr<_FuncRealToCmplx>(_ispgi_count);
    for (const auto& [espg, impedance] : _espg_to_impedance) {
        const size_t ispgi = _espg_to_ispg.at(espg);
        _ispgi_to_impedance[ispgi] = impedance;
    }
    _isei_count = 0;
    for (size_t sei=0; sei!=_sei_count; ++sei) {
        if (_espg_to_impedance.contains(_sei_to_espg[sei])) {
            ++_isei_count;
        }
    }
    _isei_to_sei = fz::SafePtr<size_t>(_isei_count);
    size_t isei = 0;
    for (size_t sei=0; sei!=_sei_count; ++sei) {
        if (_espg_to_impedance.contains(_sei_to_espg[sei])) {
            _isei_to_sei[isei] = sei;
            ++isei;
        }
    }
    _isei_to_ispgi = fz::SafePtr<size_t>(_isei_count);
    for (size_t isei=0; isei!=_isei_count; ++isei) {
        size_t sei = _isei_to_sei[isei];
        _isei_to_ispgi[isei] = _espg_to_ispg.at(_sei_to_espg[sei]);
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_organize_velocity_physical_group_data()
{
    _ispgv_to_velocity = fz::SafePtr<_FuncRealToCmplx>(_ispgv_count);
    for (const auto& [espg, volvel] : _espg_to_velocity) {
        const size_t ispgv = _espg_to_ispg.at(espg);
        _ispgv_to_velocity[ispgv] = volvel;
    }
    _vsei_count = 0;
    for (size_t sei=0; sei!=_sei_count; ++sei) {
        if (_espg_to_velocity.contains(_sei_to_espg[sei])) {
            ++_vsei_count;
        }
    }
    _vsei_to_sei = fz::SafePtr<size_t>(_vsei_count);
    size_t vsei = 0;
    for (size_t sei=0; sei!=_sei_count; ++sei) {
        if (_espg_to_velocity.contains(_sei_to_espg[sei])) {
            _vsei_to_sei[vsei] = sei;
            ++vsei;
        }
    }
    _vsei_to_ispgv = fz::SafePtr<size_t>(_vsei_count);
    for (size_t vsei=0; vsei!=_vsei_count; ++vsei) {
        size_t sei = _vsei_to_sei[vsei];
        _vsei_to_ispgv[vsei] = _espg_to_ispg.at(_sei_to_espg[sei]);
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_organize_pressure_physical_group_data()
{
    _ispgp_to_pressure = fz::SafePtr<_FuncRealToCmplx>(_ispgp_count);
    for (const auto& [espg, pressure] : _espg_to_pressure) {
        const size_t ispgp = _espg_to_ispg.at(espg);
        _ispgp_to_pressure[ispgp] = pressure;
    }
    _psei_count = 0;
    for (size_t sei=0; sei!=_sei_count; ++sei) {
        if (_espg_to_pressure.contains(_sei_to_espg[sei])) {
            ++_psei_count;
        }
    }
    _psei_to_sei = fz::SafePtr<size_t>(_psei_count);
    size_t psei = 0;
    for (size_t sei=0; sei!=_sei_count; ++sei) {
        if (_espg_to_pressure.contains(_sei_to_espg[sei])) {
            _psei_to_sei[psei] = sei;
            ++psei;
        }
    }
    _psei_to_ispgp = fz::SafePtr<size_t>(_psei_count);
    for (size_t psei=0; psei!=_psei_count; ++psei) {
        size_t sei = _psei_to_sei[psei];
        _psei_to_ispgp[psei] = _espg_to_ispg.at(_sei_to_espg[sei]);
    }
}

#if NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
    void print_dss_error(const _INTEGER_t& error_id) {
        error("oneMLK error code: {}", error_id);
    }
#endif

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_organize_physical_group_data() {
    _organize_volume_physical_group_data();
    _organize_impedance_physical_group_data();
    _organize_velocity_physical_group_data();
    _organize_pressure_physical_group_data();
}

#if NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
void define_onemkl_sparsity_pattern(
    _MKL_DSS_HANDLE_t& dss_handle,
    const fz::SafePtr<std::pair<size_t,size_t>>& nnz_rowcol_idx_pairs,
    const size_t ni_count,
    fz::SafePtr<_cmplx_t>& b_dense
) {
    // problem dimensions
    const MKL_INT node_count = ni_count;
    const MKL_INT nnz_count = nnz_rowcol_idx_pairs.size();
    
    // create a_col_idx and a_row_ptr
    fz::SafePtr<MKL_INT> a_col_idx(nnz_count);
    fz::SafePtr<MKL_INT> a_row_ptr(node_count + 1);
    MKL_INT* it_a_col_idx = a_col_idx.begin();
    MKL_INT* it_a_row_ptr = a_row_ptr.begin();
    const std::pair<size_t,size_t>* it_nnz_rowcol_idx_pairs =
        nnz_rowcol_idx_pairs.begin();
    size_t current_row = std::numeric_limits<size_t>::max();
    for (MKL_INT i=0; i!=nnz_count; ++i) {
        *it_a_col_idx = it_nnz_rowcol_idx_pairs->second;
        ++it_a_col_idx;
        if (it_nnz_rowcol_idx_pairs->first != current_row) {
            current_row = it_nnz_rowcol_idx_pairs->first;
            *it_a_row_ptr = i; 
            ++it_a_row_ptr;
        }
        ++it_nnz_rowcol_idx_pairs;
    }
    *it_a_row_ptr = nnz_count;
    ++it_a_row_ptr;
    
    constexpr MKL_INT options = NUMAV_MKL_OPTIONS;
    
    // error code
    _INTEGER_t error_id;
    
    // initialize the solver
    error_id = dss_create(dss_handle, options);
    if (error_id != MKL_DSS_SUCCESS) { print_dss_error(error_id); }
    
    // define the non-zero structure of the matrix
    const MKL_INT symmetry_type = MKL_DSS_SYMMETRIC_COMPLEX;
    error_id = dss_define_structure(
        dss_handle, symmetry_type, a_row_ptr.data(),
        node_count, node_count, a_col_idx.data(), nnz_count
    );
    if (error_id != MKL_DSS_SUCCESS) { print_dss_error(error_id); }
    
    // reorder the matrix
    error_id = dss_reorder(dss_handle, options, 0);
    if (error_id != MKL_DSS_SUCCESS) { print_dss_error(error_id); }
    a_row_ptr.free();
    a_col_idx.free();

    // allocate the null dense b vector
    b_dense = fz::SafePtr<_cmplx_t>(ni_count);
    b_dense.fill(0.0);
}
#endif

template<typename T>
bool compare_pair(const std::pair<T,T>& a, const std::pair<T,T>& b) {
    #if NUMAV_GLOBAL_MATRIX_STORAGE_ORDER == NUMAV_UPPER_COL_MAJOR
        if (a.second != b.second) {
            return a.second < b.second;
        }
        else {
            return a.first < b.first;
        }
    #elif NUMAV_GLOBAL_MATRIX_STORAGE_ORDER == NUMAV_UPPER_ROW_MAJOR
        if (a.first != b.first) {
            return a.first < b.first;
        }
        else {
            return a.second < b.second;
        }
    #else
        static_assert(false, "Invalid GLOBAL_MATRIX_STORAGE_ORDER.");
    #endif
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_allocate_a_and_b()
{
    // allocate a
    constexpr std::array<
        std::array<size_t,2>, COMB_REP_SIZE<NODES_IN_VOL_ELEM<O>,2>
    > COMBS_VOL = COMBINATION_REP<NODES_IN_VOL_ELEM<O>>;
    std::unordered_set<std::pair<size_t,size_t>> existing_pairs;
    for (size_t vei=0; vei!=_vei_count; ++vei) {
        for (const auto& c : COMBS_VOL) {
            existing_pairs.insert( 
                make_ascending_pair(
                    _vei_to_ni[vei][c[0]], _vei_to_ni[vei][c[1]]
                )
            );
        }
    }
    _nnz_rowcol_idx_pairs = fz::SafePtr<std::pair<size_t,size_t>>(
        existing_pairs.size()
    );
    std::copy(
        existing_pairs.begin(),
        existing_pairs.end(),
        _nnz_rowcol_idx_pairs.begin()
    );
    std::sort(
        _nnz_rowcol_idx_pairs.begin(),
        _nnz_rowcol_idx_pairs.end(),
        compare_pair<size_t>
    );
    _a_vals = fz::SafePtr<_cmplx_t>(_nnz_rowcol_idx_pairs.size());

    // allocate b
    std::unordered_set<size_t> existing_source_nodes;
    for (size_t vsei=0; vsei!=_vsei_count; ++vsei) {
        const size_t sei = _vsei_to_sei[vsei];
        for (size_t eni=0; eni!=NODES_IN_SFC_ELEM<O>; ++eni) {
            const size_t ni = _sei_to_ni[sei][eni];
            existing_source_nodes.insert(ni);
        }
    }
    for (size_t pvni=0; pvni!=_pvni_count; ++pvni) {
        existing_source_nodes.insert(std::get<size_t>(_point_volvel[pvni]));
    }
    for (size_t psei=0; psei!=_psei_count; ++psei) {
        const size_t sei = _psei_to_sei[psei];
        for (size_t eni=0; eni!=NODES_IN_SFC_ELEM<O>; ++eni) {
            const size_t ni = _sei_to_ni[sei][eni];
            existing_source_nodes.insert(ni);
        }
    }
    for (size_t ppni=0; ppni!=_ppni_count; ++ppni) {
        existing_source_nodes.insert(std::get<size_t>(_point_pressure[ppni]));
    }
    _b_row_idx = fz::SafePtr<size_t>(existing_source_nodes.size());
    std::copy(
        existing_source_nodes.begin(),
        existing_source_nodes.end(),
        _b_row_idx.begin()
    );
    std::sort(_b_row_idx.begin(), _b_row_idx.end());
    _b_vals = fz::SafePtr<_cmplx_t>(_b_row_idx.size());
}

template<ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_assemble_fi_part_for_vol_elements()
{
    constexpr std::array<
        std::array<size_t,2>, COMB_REP_SIZE<NODES_IN_VOL_ELEM<O>,2>
    > COMBS_VOL = COMBINATION_REP<NODES_IN_VOL_ELEM<O>>;

    // count the fipi for each ivpg
    fz::SafePtr<std::unordered_map<
        std::pair<size_t,size_t>, size_t
    >> ivpg_to_map_to_fipi(_ivpg_count);
    for (size_t vei=0; vei!=_vei_count; ++vei)
    {
        const size_t ivpg = _vei_to_ivpg[vei];
        for (const auto& c : COMBS_VOL) {
            const std::pair<size_t,size_t> pair = make_ascending_pair(
                _vei_to_ni[vei][c[0]], _vei_to_ni[vei][c[1]]
            );
            if (!ivpg_to_map_to_fipi[ivpg].contains(pair)) {
                const size_t fipi = ivpg_to_map_to_fipi[ivpg].size();
                ivpg_to_map_to_fipi[ivpg].insert({pair, fipi});
            }
        }
    }

    // allocate memory in the safe pointers
    _ivpg_to_stif_fi_part = fz::SafePtr<fz::SafePtr<double>>(_ivpg_count);
    _ivpg_to_mass_fi_part = fz::SafePtr<fz::SafePtr<double>>(_ivpg_count);
    _ivpg_to_ptr_in_a = fz::SafePtr<fz::SafePtr<_cmplx_t*>>(_ivpg_count);
    for (size_t ivpg=0; ivpg!=_ivpg_count; ++ivpg)
    {
        const size_t size = ivpg_to_map_to_fipi[ivpg].size();
        _ivpg_to_stif_fi_part[ivpg] = fz::SafePtr<double>(size);
        _ivpg_to_stif_fi_part[ivpg].fill(0.0);
        _ivpg_to_mass_fi_part[ivpg] = fz::SafePtr<double>(size);
        _ivpg_to_mass_fi_part[ivpg].fill(0.0);
        _ivpg_to_ptr_in_a[ivpg] = fz::SafePtr<_cmplx_t*>(size);
    }

    // assemble elementary stiffness and mass matrices
    std::array<size_t, COMBS_VOL.size()> fipi_vol;
    for (size_t vei=0; vei!=_vei_count; ++vei)
    {
        const size_t ivpg = _vei_to_ivpg[vei];
        
        // Create _ivpg_to_ptr_in_a
        for (size_t nci=0; nci!=COMBS_VOL.size(); ++nci)
        {
            const std::pair<size_t,size_t> pair = make_ascending_pair(
                _vei_to_ni[vei][COMBS_VOL[nci][0]],
                _vei_to_ni[vei][COMBS_VOL[nci][1]]
            );
            fipi_vol[nci] = ivpg_to_map_to_fipi[ivpg].at(pair);
            
            const std::pair<size_t,size_t>* const pair_ptr = std::lower_bound(
                _nnz_rowcol_idx_pairs.begin(),
                _nnz_rowcol_idx_pairs.end(),
                pair,
                compare_pair<size_t>
            );
            const ptrdiff_t ptrdiff = pair_ptr - _nnz_rowcol_idx_pairs.begin();
            _ivpg_to_ptr_in_a[ivpg][fipi_vol[nci]] = _a_vals.begin() + ptrdiff;
        }
        
        // coordinates matrix
        Eigen::Matrix<double,NODES_IN_VOL_ELEM<O>,DIM> coords_matrix;
        for (size_t ni=0; ni!=NODES_IN_VOL_ELEM<O>; ++ni) {
            const size_t node_idx = _vei_to_ni[vei][ni];
            for (size_t di=0; di!=DIM; ++di) {
                coords_matrix(ni,di) = _node_coords[node_idx][di];
            }
        }

        // stiffness matrix
        constexpr std::array<std::array<double,DIM>,NGP_STIF<O>>
            GAUSS_POINTS_STIF = GAUSS_POINTS_VOL<NGP_STIF<O>>;
        for (size_t gpi=0; gpi!=NGP_STIF<O>; ++gpi)
        {
            const Eigen::Matrix<double,DIM,NODES_IN_VOL_ELEM<O>> nabla_n =
                shape_func_vol_gradient<O>(
                    GAUSS_POINTS_STIF[gpi][0],
                    GAUSS_POINTS_STIF[gpi][1],
                    GAUSS_POINTS_STIF[gpi][2]
                ); // todo: try putting constexpr here

            const Eigen::Matrix<double,DIM,DIM> jac_matrix = 
                nabla_n * coords_matrix;

            const Eigen::Matrix<double,DIM,DIM> inv_jac = jac_matrix.inverse();
            
            const Eigen::Matrix<double,DIM,NODES_IN_VOL_ELEM<O>> b_matrix =
                inv_jac * nabla_n;
            
            const
            Eigen::Matrix<double,NODES_IN_VOL_ELEM<O>,NODES_IN_VOL_ELEM<O>>
                btb = b_matrix.transpose() * b_matrix;
            
            const double det_jac = jac_matrix.determinant();
            
            // todo: multiply detj and w without creating another eigen matrix
            const
            Eigen::Matrix<double,NODES_IN_VOL_ELEM<O>,NODES_IN_VOL_ELEM<O>>
                btb_detj_w = 
                    btb * det_jac * GAUSS_WEIGHTS_VOL<NGP_STIF<O>>[gpi];

            for (size_t nci=0; nci!=COMBS_VOL.size(); ++nci) {
                _ivpg_to_stif_fi_part[ivpg][fipi_vol[nci]] += btb_detj_w(
                    COMBS_VOL[nci][0], COMBS_VOL[nci][1]
                );
            }
        }

        // mass matrix
        constexpr std::array<std::array<double,DIM>,NGP_MASS<O>>
            GAUSS_POINTS_MASS = GAUSS_POINTS_VOL<NGP_MASS<O>>;
        for (size_t gpi=0; gpi!=NGP_MASS<O>; ++gpi)
        {
            const Eigen::Matrix<double,DIM,NODES_IN_VOL_ELEM<O>> nabla_n =
                shape_func_vol_gradient<O>(
                    GAUSS_POINTS_MASS[gpi][0],
                    GAUSS_POINTS_MASS[gpi][1],
                    GAUSS_POINTS_MASS[gpi][2]
                ); // todo: try putting constexpr here

            const Eigen::Matrix<double,DIM,DIM> jac_matrix = 
                nabla_n * coords_matrix;
            
            const double det_jac = jac_matrix.determinant();

            const Eigen::Matrix<double,NODES_IN_VOL_ELEM<O>,1> n =
                shape_func_vol<O>(
                    GAUSS_POINTS_MASS[gpi][0],
                    GAUSS_POINTS_MASS[gpi][1],
                    GAUSS_POINTS_MASS[gpi][2]
                );

            const
            Eigen::Matrix<double,NODES_IN_VOL_ELEM<O>,NODES_IN_VOL_ELEM<O>>
                nnt = n * n.transpose();
            
            // todo: multiply detj and w without creating another eigen matrix
            const
            Eigen::Matrix<double,NODES_IN_VOL_ELEM<O>,NODES_IN_VOL_ELEM<O>>
                nnt_detj_w =
                    nnt * det_jac * GAUSS_WEIGHTS_VOL<NGP_MASS<O>>[gpi];
            
            for (size_t nci=0; nci!=COMBS_VOL.size(); ++nci) {
                _ivpg_to_mass_fi_part[ivpg][fipi_vol[nci]] += nnt_detj_w(
                    COMBS_VOL[nci][0], COMBS_VOL[nci][1]
                );
            }
        }
    }
    ivpg_to_map_to_fipi.free();
}

template<ElementOrder O>
std::array<Eigen::Vector2d, NODES_IN_SFC_ELEM<O>> project_triangle_to_2d(
    const std::array<Eigen::Vector3d, NODES_IN_SFC_ELEM<O>>& vertices_3d
) {
    std::array<Eigen::Vector3d, NODES_IN_SFC_ELEM<O>> point_minus_origin;
    for (size_t ni = 0; ni!=NODES_IN_SFC_ELEM<O>; ++ni) {
        point_minus_origin[ni] = vertices_3d[ni] - vertices_3d[0];
    }
    const Eigen::Vector3d& u = point_minus_origin[1];
    const Eigen::Vector3d& v = point_minus_origin[2];
    const Eigen::Vector3d n = u.cross(v);
    const Eigen::Vector3d x = u / u.norm();
    Eigen::Vector3d y = n.cross(u);
    y /= y.norm();

    std::array<Eigen::Vector2d, NODES_IN_SFC_ELEM<O>> vertices_2d;
    for (size_t ni = 0; ni!=NODES_IN_SFC_ELEM<O>; ++ni) {
        vertices_2d[ni] = Eigen::Vector2d(
            point_minus_origin[ni].dot(x), point_minus_origin[ni].dot(y)
        );
    }
    return vertices_2d;
}

template<ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_assemble_fi_part_for_sfc_impedance()
{
    constexpr std::array<
        std::array<size_t,2>, COMB_REP_SIZE<NODES_IN_SFC_ELEM<O>,2>
    > COMBS_SFC = COMBINATION_REP<NODES_IN_SFC_ELEM<O>>;

    // count the fipi for each ispgi
    fz::SafePtr<std::unordered_map<
        std::pair<size_t,size_t>, size_t
    >> ispgi_to_map_to_fipi(_ispgi_count);
    for (size_t isei=0; isei!=_isei_count; ++isei)
    {
        const size_t ispgi = _isei_to_ispgi[isei];
        const size_t sei = _isei_to_sei[isei];
        for (const auto& c : COMBS_SFC) {
            const std::pair<size_t,size_t> pair = make_ascending_pair(
                _sei_to_ni[sei][c[0]], _sei_to_ni[sei][c[1]]
            );
            if (!ispgi_to_map_to_fipi[ispgi].contains(pair)) {
                const size_t fipi = ispgi_to_map_to_fipi[ispgi].size();
                ispgi_to_map_to_fipi[ispgi].insert({pair, fipi});
            }
        }
    }

    // allocate memory in the safe pointers
    _ispgi_to_damp_fi_part = fz::SafePtr<fz::SafePtr<double>>(_ispgi_count);
    _ispgi_to_ptr_in_a = fz::SafePtr<fz::SafePtr<_cmplx_t*>>(_ispgi_count);
    for (size_t ispgi=0; ispgi!=_ispgi_count; ++ispgi)
    {
        const size_t size = ispgi_to_map_to_fipi[ispgi].size();
        _ispgi_to_damp_fi_part[ispgi] = fz::SafePtr<double>(size);
        _ispgi_to_damp_fi_part[ispgi].fill(0.0);
        _ispgi_to_ptr_in_a[ispgi] = fz::SafePtr<_cmplx_t*>(size);
    }

    // assemble the elementary damping matrices
    std::array<size_t, COMBS_SFC.size()> fipi_damp;
    for (size_t isei=0; isei!=_isei_count; ++isei)
    {
        const size_t ispgi = _isei_to_ispgi[isei];
        const size_t sei = _isei_to_sei[isei];
        
        // Create _ispgi_to_ptr_in_a
        for (size_t nci=0; nci!=COMBS_SFC.size(); ++nci)
        {
            const std::pair<size_t,size_t> pair = make_ascending_pair(
                _sei_to_ni[sei][COMBS_SFC[nci][0]],
                _sei_to_ni[sei][COMBS_SFC[nci][1]]
            );
            fipi_damp[nci] = ispgi_to_map_to_fipi[ispgi].at(pair);
            
            const std::pair<size_t,size_t>* const pair_ptr = std::lower_bound(
                _nnz_rowcol_idx_pairs.begin(),
                _nnz_rowcol_idx_pairs.end(),
                pair,
                compare_pair<size_t>
            );
            const ptrdiff_t ptrdiff = pair_ptr - _nnz_rowcol_idx_pairs.begin();
            _ispgi_to_ptr_in_a[ispgi][fipi_damp[nci]] =
                _a_vals.begin() + ptrdiff;
        }

        #if NUMAV_TRIANGLE_INTEGRATION_METHOD == NUMAV_JACOBIAN_DETERMINANT
            // coordinates matrix
            std::array<Eigen::Vector3d, NODES_IN_SFC_ELEM<O>> triangle_3d;
            for (size_t eni=0; eni!=NODES_IN_SFC_ELEM<O>; ++eni) {
                const size_t node_idx = _sei_to_ni[sei][eni];
                triangle_3d[eni] = Eigen::Vector3d(
                    _node_coords[node_idx][0],
                    _node_coords[node_idx][1],
                    _node_coords[node_idx][2]
                );
            }
            std::array<Eigen::Vector2d, NODES_IN_SFC_ELEM<O>> triangle_2d =
                project_triangle_to_2d<O>(triangle_3d);
            Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,2> coords_matrix;
            for (size_t eni=0; eni!=NODES_IN_SFC_ELEM<O>; ++eni) {
                coords_matrix(eni,0) = triangle_2d[eni](0);
                coords_matrix(eni,1) = triangle_2d[eni](1);
            }
        #elif NUMAV_TRIANGLE_INTEGRATION_METHOD == NUMAV_TRIANGLE_AREA
            // Calculate triangle area
            std::array<std::array<double,DIM>,3> triangle_coords;
            for (size_t ni=0; ni!=3; ++ni) {
                const size_t node_idx = _sei_to_ni[sei][ni];
                triangle_coords[ni] = std::array<double,DIM>({
                    _node_coords[node_idx][0],
                    _node_coords[node_idx][1],
                    _node_coords[node_idx][2]
                });
            }
            const double det_jac = 2*get_triangle_area(triangle_coords);
        #else
            static_assert(
                false, "Invalid NUMAV_TRIANGLE_INTEGRATION_METHOD."
            );
        #endif

        // damping matrix
        constexpr std::array<std::array<double,2>,NGP_DAMP<O>>
            GAUSS_POINTS_DAMP = GAUSS_POINTS_SFC<NGP_DAMP<O>>;
        for (size_t gpi=0; gpi!=NGP_DAMP<O>; ++gpi)
        {
            #if NUMAV_TRIANGLE_INTEGRATION_METHOD == NUMAV_JACOBIAN_DETERMINANT
                const Eigen::Matrix<double,2,NODES_IN_SFC_ELEM<O>> nabla_n =
                    shape_func_sfc_gradient<O>(
                        GAUSS_POINTS_DAMP[gpi][0], GAUSS_POINTS_DAMP[gpi][1]
                    ); // todo: try putting constexpr here
                const Eigen::Matrix<double,2,2> jac_matrix =
                    nabla_n * coords_matrix;
                const double det_jac = jac_matrix.determinant();
            #endif

            const Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,1> n =
                shape_func_sfc<O>(
                    GAUSS_POINTS_DAMP[gpi][0], GAUSS_POINTS_DAMP[gpi][1]
                );
            
            const
            Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,NODES_IN_SFC_ELEM<O>>
                nnt = n * n.transpose();

            // todo: multiply detj and w without creating another eigen matrix
            const
            Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,NODES_IN_SFC_ELEM<O>>
                nnt_detj_w =
                    nnt * det_jac * GAUSS_WEIGHTS_SFC<NGP_DAMP<O>>[gpi];
            
            for (size_t nci=0; nci!=COMBS_SFC.size(); ++nci) {
                _ispgi_to_damp_fi_part[ispgi][fipi_damp[nci]] += nnt_detj_w(
                    COMBS_SFC[nci][0], COMBS_SFC[nci][1]
                );
            }
        }
    }
    ispgi_to_map_to_fipi.free();
}

template<ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_assemble_fi_part_for_point_velocity()
{
    _pvni_to_forc_fi_part = fz::SafePtr<_FuncRealToCmplx>(_pvni_count);
    _pvni_to_ptr_in_b = fz::SafePtr<_cmplx_t*>(_pvni_count);
    for (size_t pvni=0; pvni!=_pvni_count; ++pvni)
    {
        _pvni_to_forc_fi_part[pvni] = 
            std::get<_FuncRealToCmplx>(_point_volvel[pvni]);

        const size_t ni = std::get<size_t>(_point_volvel[pvni]);
        const size_t* const ni_ptr = std::lower_bound(
            _b_row_idx.begin(), _b_row_idx.end(), ni
        );
        const ptrdiff_t ptrdiff = ni_ptr - _b_row_idx.begin();
        _pvni_to_ptr_in_b[pvni] = _b_vals.begin() + ptrdiff;
    }
}

template<ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_assemble_fi_part_for_sfc_velocity()
{
    // count the fipi for each ispgv
    fz::SafePtr<std::unordered_map<size_t,size_t>> ispgv_to_map_to_fipi(
        _ispgv_count
    );
    for (size_t vsei=0; vsei!=_vsei_count; ++vsei) {
        const size_t ispgv = _vsei_to_ispgv[vsei];
        const size_t sei = _vsei_to_sei[vsei];
        for (size_t eni=0; eni!=NODES_IN_SFC_ELEM<O>; ++eni) {
            const size_t ni = _sei_to_ni[sei][eni];
            if (!ispgv_to_map_to_fipi[ispgv].contains(ni)) {
                const size_t fipi = ispgv_to_map_to_fipi[ispgv].size();
                ispgv_to_map_to_fipi[ispgv].insert({ni, fipi});
            }
        }
    }

    // allocate memory in the safe pointers
    _ispgv_to_forc_fi_part = fz::SafePtr<fz::SafePtr<double>>(_ispgv_count);
    _ispgv_to_ptr_in_b = fz::SafePtr<fz::SafePtr<_cmplx_t*>>(_ispgv_count);
    for (size_t ispgv=0; ispgv!=_ispgv_count; ++ispgv)
    {
        const size_t size = ispgv_to_map_to_fipi[ispgv].size();
        _ispgv_to_forc_fi_part[ispgv] = fz::SafePtr<double>(size);
        _ispgv_to_forc_fi_part[ispgv].fill(0.0);
        _ispgv_to_ptr_in_b[ispgv] = fz::SafePtr<_cmplx_t*>(size);
    }

    // assemble elementary force vectors
    std::array<size_t, NODES_IN_SFC_ELEM<O>> fipi_forc;
    for (size_t vsei=0; vsei!=_vsei_count; ++vsei)
    {
        const size_t ispgv = _vsei_to_ispgv[vsei];
        const size_t sei = _vsei_to_sei[vsei];

        // Create _ispgv_to_ptr_in_b
        for (size_t eni=0; eni!=NODES_IN_SFC_ELEM<O>; ++eni)
        {
            const size_t ni = _sei_to_ni[sei][eni];
            fipi_forc[eni] = ispgv_to_map_to_fipi[ispgv].at(ni);
            
            const size_t* const ni_ptr = std::lower_bound(
                _b_row_idx.begin(), _b_row_idx.end(), ni
            );
            const ptrdiff_t ptrdiff = ni_ptr - _b_row_idx.begin();
            _ispgv_to_ptr_in_b[ispgv][fipi_forc[eni]] =
                _b_vals.begin() + ptrdiff;
        }

        // 2D coordinates matrix
        std::array<Eigen::Vector3d, NODES_IN_SFC_ELEM<O>> triangle_3d;
        for (size_t eni=0; eni!=NODES_IN_SFC_ELEM<O>; ++eni) {
            const size_t node_idx = _sei_to_ni[sei][eni];
            triangle_3d[eni] = Eigen::Vector3d(
                _node_coords[node_idx][0],
                _node_coords[node_idx][1],
                _node_coords[node_idx][2]
            );
        }
        std::array<Eigen::Vector2d, NODES_IN_SFC_ELEM<O>> triangle_2d =
            project_triangle_to_2d<O>(triangle_3d);
        Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,2> coords_matrix;
        for (size_t eni=0; eni!=NODES_IN_SFC_ELEM<O>; ++eni) {
            coords_matrix(eni,0) = triangle_2d[eni](0);
            coords_matrix(eni,1) = triangle_2d[eni](1);
        }

        // elementary force vector
        constexpr std::array<std::array<double,2>,NGP_FORC<O>>
            GAUSS_POINTS_FORC = GAUSS_POINTS_SFC<NGP_FORC<O>>;
        for (size_t gpi=0; gpi!=NGP_FORC<O>; ++gpi)
        {
            #if NUMAV_TRIANGLE_INTEGRATION_METHOD == NUMAV_JACOBIAN_DETERMINANT
                const Eigen::Matrix<double,2,NODES_IN_SFC_ELEM<O>> nabla_n =
                    shape_func_sfc_gradient<O>(
                        GAUSS_POINTS_FORC[gpi][0], GAUSS_POINTS_FORC[gpi][1]
                    ); // todo: try putting constexpr here
                const Eigen::Matrix<double,2,2> jac_matrix =
                    nabla_n * coords_matrix;
                const double det_jac = jac_matrix.determinant();
            #endif

            const Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,1> n =
                shape_func_sfc<O>(
                    GAUSS_POINTS_FORC[gpi][0], GAUSS_POINTS_FORC[gpi][1]
                );

            // todo: multiply detj and w without creating another eigen matrix
            const Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,1> n_detj_w =
                n * det_jac * GAUSS_WEIGHTS_SFC<NGP_FORC<O>>[gpi];
            
            for (size_t eni=0; eni!=NODES_IN_SFC_ELEM<O>; ++eni) {
                _ispgv_to_forc_fi_part[ispgv][fipi_forc[eni]] += n_detj_w(eni);
            }
        }
    }
    ispgv_to_map_to_fipi.free();
}

template<typename T>
std::unordered_map<std::vector<size_t>, std::vector<T>> find_set_intersections(
    const fz::SafePtr<fz::SafePtr<T>>& sets
) {
    // map each element to the indices of sets that contain it
    std::unordered_map<T, std::vector<size_t>> element_to_sets;
    for (size_t set_index=0; set_index!=sets.size(); ++set_index) {
        for (const auto& element : sets[set_index]) {
            element_to_sets[element].push_back(set_index);
        }
    }
    // group elements by which sets contain them
    std::unordered_map<std::vector<size_t>, std::vector<T>> intersections;
    for (const auto& [element, set_indices] : element_to_sets) {
        std::vector<size_t> sorted_indices = set_indices;
        std::sort(sorted_indices.begin(), sorted_indices.end());
        intersections[sorted_indices].push_back(element);
    }
    return intersections;
}

template<ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_assemble_fi_part_for_pressure()
{
    // create the mathmatical sets of nodes for each pressure value assigned
    fz::SafePtr<fz::SafePtr<size_t>> sets(_ppni_count + _ispgp_count);
    for (size_t ppni=0; ppni!=_ppni_count; ++ppni) {
        sets[ppni] = fz::SafePtr<size_t>(1);
        sets[ppni][0] = std::get<size_t>(_point_pressure[ppni]);
    }
    for (size_t ispgp=0; ispgp!=_ispgp_count; ++ispgp) {
        std::set<size_t> unique_nodes;
        for (size_t psei=0; psei!=_psei_count; ++psei) {
            if (_psei_to_ispgp[psei] == ispgp){
                const size_t sei = _psei_to_sei[psei];
                for (size_t eni=0; eni!=NODES_IN_SFC_ELEM<O>; ++eni) {
                    const size_t ni = _sei_to_ni[sei][eni];
                    unique_nodes.insert(ni);
                }
            }
        }
        sets[_ppni_count+ispgp] = fz::SafePtr<size_t>(unique_nodes.size());
        size_t i = 0;
        for (const auto& ni : unique_nodes) {
            sets[_ppni_count+ispgp][i] = ni;
            ++i;
        }
    }
    std::unordered_map<std::vector<size_t>,std::vector<size_t>>
        intersections = find_set_intersections(sets);
    for (auto& set : sets) {
        set.free();
    }
    sets.free();

    // calculate the average pressure between intersected elements in the sets
    _pvi_count = intersections.size();
    _pvi_to_pressure = fz::SafePtr<_FuncRealToCmplx>(_pvi_count);
    size_t pvi = 0;
    for (auto& [set_indexes, ni_vector] : intersections) {
        auto average = [this,set_indexes](const double& freq) {
            _cmplx_t sum = 0;
            for (const auto& set_index : set_indexes) {
                if (set_index < _ppni_count) {
                    const size_t& ppni = set_index;
                    sum += std::get<_FuncRealToCmplx>(
                        _point_pressure[ppni]
                    )(freq);
                }
                else {
                    const size_t ispgp = set_index - _ppni_count;
                    sum += (_ispgp_to_pressure[ispgp])(freq);
                }
            }
            return sum / static_cast<double>(set_indexes.size());
        };
        _pvi_to_pressure[pvi] = average;
        ++pvi;
    }

    // define _pvi_to_ptr_in_a and _pvi_to_ptr_in_b
    _pvi_to_ptr_in_a = fz::SafePtr<fz::SafePtr<_cmplx_t*>>(_pvi_count);
    _pvi_to_ptr_in_b = fz::SafePtr<fz::SafePtr<_cmplx_t*>>(_pvi_count);
    pvi = 0;
    for (auto& [set_indexes, ni_vector] : intersections) {
        _pvi_to_ptr_in_a[pvi] = fz::SafePtr<_cmplx_t*>(ni_vector.size());
        _pvi_to_ptr_in_b[pvi] = fz::SafePtr<_cmplx_t*>(ni_vector.size());
        size_t fipi = 0;
        for (const auto& ni : ni_vector) {
            // _pvi_to_ptr_in_a
            const std::pair<size_t,size_t> pair = std::make_pair(ni, ni);
            const std::pair<size_t,size_t>* const pair_ptr =
                std::lower_bound(
                    _nnz_rowcol_idx_pairs.begin(),
                    _nnz_rowcol_idx_pairs.end(),
                    pair,
                    compare_pair<size_t>
                );
            const ptrdiff_t ptrdiff_a = 
                pair_ptr - _nnz_rowcol_idx_pairs.begin();
            _pvi_to_ptr_in_a[pvi][fipi] = _a_vals.begin() + ptrdiff_a;
            
            // _pvi_to_ptr_in_b
            const size_t* const ni_ptr = std::lower_bound(
                _b_row_idx.begin(), _b_row_idx.end(), ni
            );
            const ptrdiff_t ptrdiff_b = ni_ptr - _b_row_idx.begin();
            _pvi_to_ptr_in_b[pvi][fipi] = _b_vals.begin() + ptrdiff_b;
            
            ++fipi;
        }
        ++pvi;
    }
}

template<ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_assemble_freq_independent_parts()
{   
    #if NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
        define_onemkl_sparsity_pattern(
            _dss_handle, _nnz_rowcol_idx_pairs, _ni_count, _b_dense
        );
    #endif
    _assemble_fi_part_for_vol_elements();
    _assemble_fi_part_for_sfc_impedance();
    _assemble_fi_part_for_point_velocity();
    _assemble_fi_part_for_sfc_velocity();
    _assemble_fi_part_for_pressure();
}

#if NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
void solve_using_onemkl(
    _MKL_DSS_HANDLE_t& dss_handle,
    const fz::SafePtr<_cmplx_t>& a_vals,
    const fz::SafePtr<_cmplx_t>& b_vals,
    const fz::SafePtr<size_t>& b_row_idx,
    fz::SafePtr<_cmplx_t>& b_dense,
    _cmplx_t* const x_out
) {
    // error code
    _INTEGER_t error_id;
    
    // factor the matrix
    constexpr MKL_INT positive_definiteness = MKL_DSS_INDEFINITE;
    error_id = dss_factor_complex(
        dss_handle, positive_definiteness,
        reinterpret_cast<const double*>(a_vals.data())
    );
    if (error_id != MKL_DSS_SUCCESS) { print_dss_error(error_id); }

    // dense b vector
    for (size_t i=0; i!=b_vals.size(); ++i) {
        b_dense[b_row_idx[i]] = b_vals[i];
    }
    
    // solve
    constexpr MKL_INT options = NUMAV_MKL_OPTIONS;
    constexpr MKL_INT num_of_b = 1;
    error_id = dss_solve_complex(
        dss_handle, options, reinterpret_cast<const double*>(b_dense.data()),
        num_of_b, reinterpret_cast<double*>(x_out)
    );
    if (error_id != MKL_DSS_SUCCESS) { print_dss_error(error_id); }
}
#endif

void solve_using_eigen(
    const fz::SafePtr<_cmplx_t>& a_vals,
    const fz::SafePtr<std::pair<size_t,size_t>>& nnz_rowcol_idx_pairs,
    const fz::SafePtr<_cmplx_t>& b_vals,
    const fz::SafePtr<size_t>& b_row_idx,
    const size_t& node_count,
    _cmplx_t* const x_out
) {
    using Triplet = typename Eigen::Triplet<_cmplx_t>;

    // a matrix
    fz::SafePtr<Triplet> triplets_a(2*a_vals.size() - node_count);
    Triplet* it_triplets_a = triplets_a.begin();
    for (size_t j=0; j!=a_vals.size(); ++j) {
        *it_triplets_a = Triplet(
            nnz_rowcol_idx_pairs[j].first,
            nnz_rowcol_idx_pairs[j].second,
            a_vals[j]
        );
        ++it_triplets_a;

        // lower part of a
        if (nnz_rowcol_idx_pairs[j].first != 
            nnz_rowcol_idx_pairs[j].second
        ) {
            *it_triplets_a = Triplet(
                nnz_rowcol_idx_pairs[j].second,
                nnz_rowcol_idx_pairs[j].first,
                a_vals[j]
            );
            ++it_triplets_a;
        }
    }
    Eigen::SparseMatrix<_cmplx_t> a(node_count, node_count);
    a.setFromTriplets(triplets_a.begin(), triplets_a.end());
    triplets_a.free();

    // b vector
    fz::SafePtr<Triplet> triplets_b(node_count);
    Triplet* it_triplets_b = triplets_b.begin();
    for (size_t j=0; j!=b_vals.size(); ++j) {
        *it_triplets_b = Triplet(b_row_idx[j], 0, b_vals[j]);
        ++it_triplets_b;
    }
    Eigen::SparseMatrix<_cmplx_t> b(node_count, 1);
    b.setFromTriplets(triplets_b.begin(), triplets_b.end());
    triplets_b.free();

    Eigen::SparseLU<Eigen::SparseMatrix<_cmplx_t>> solver;
    solver.analyzePattern(a);
    solver.factorize(a);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Factorization failed. Matrix may be singular.\n";
    }
    const Eigen::VectorXcd x = solver.solve(b);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solving failed.\n";
    }
    for (size_t j=0; j!=node_count; ++j) {
        x_out[j] = x(j);
    }
}

template<ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_solve()
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

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::run()
{
    _check_if_it_can_run();
    log::print_opening();
    log::print_opening_ac_fem_freq_d3();
    _define_freq_vector();
    _organize_physical_group_data();
    _allocate_a_and_b();
    _assemble_freq_independent_parts();
    _solve();
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::export_result(
    const char* const path_to_result
) {
    if (!_did_run) {
        error("The Simulation needs to run before the result is exported."
            " Call the run method first.");
    }
    write_matrix(_cmplx_pressure_amp, path_to_result);
}

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
