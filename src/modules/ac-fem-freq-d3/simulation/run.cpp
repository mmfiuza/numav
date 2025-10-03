// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav/numav.hpp"

#include <tuple>
#include <cmath>
#include <charconv>
#include <fstream>
#include <cassert>
#include <numbers>
#include <algorithm>
#include <limits>
#include <chrono>
#include <iomanip>

#include "common/log.hpp"
#include "common/hash-functions.hpp"
#include "common/maths.hpp"
#include "common/utils.hpp"

#include "mkl_dss.h"
#include "mkl_types.h"

namespace numav {

static constexpr size_t DIM = DIM_COUNT<Dimension::D3>;
static constexpr double AREA_REF_TRIG = 1.0 / 2.0;
static constexpr double VOLUME_REF_TET = 1.0 / 6.0;

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
    if constexpr (N == 15) {
        constexpr double a = 1.0 / 4.0;
        constexpr double b = ( 7.0 +     std::sqrt(15.0)) / 34.0;
        constexpr double c = ( 7.0 -     std::sqrt(15.0)) / 34.0;
        constexpr double d = (13.0 - 3.0*std::sqrt(15.0)) / 34.0;
        constexpr double e = (13.0 + 3.0*std::sqrt(15.0)) / 34.0;
        constexpr double f = ( 5.0 -     std::sqrt(15.0)) / 20.0;
        constexpr double g = ( 5.0 +     std::sqrt(15.0)) / 20.0;
        return std::array<std::array<double,DIM>,N>{{
            {a,a,a}, {b,b,b}, {b,b,d}, {b,d,b}, {d,b,b}, 
            {c,c,c}, {c,c,e}, {c,e,c}, {e,c,c}, {f,f,g},
            {f,g,f}, {g,f,f}, {f,g,g}, {g,f,g}, {g,g,f}
        }};
    }
}();

template<size_t N>
constexpr std::array<double,N> GAUSS_WEIGHTS_SFC = [] {
    if constexpr (N == 1) {
        constexpr double a = AREA_REF_TRIG;
        return std::array<double,N>({a});
    }
    if constexpr (N == 3) {
        constexpr double a = AREA_REF_TRIG * 1.0 / 3.0;
        return std::array<double,N>({a,a,a});
    }
}();

template<size_t N>
constexpr std::array<double,N> GAUSS_WEIGHTS_VOL = [] {
    if constexpr (N == 1) {
        constexpr double a = VOLUME_REF_TET;
        return std::array<double,N>({a});
    }
    if constexpr (N == 4) {
        constexpr double a = VOLUME_REF_TET * 1.0 / 4.0;
        return std::array<double,N>({a,a,a,a});
    }
    if constexpr (N == 15) {
        constexpr double a =
            VOLUME_REF_TET * 16.0 / 135.0;
        constexpr double b =
            VOLUME_REF_TET * (2665.0 - 14.0*std::sqrt(15.0)) / 37800.0;
        constexpr double c =
            VOLUME_REF_TET * (2665.0 + 14.0*std::sqrt(15.0)) / 37800.0;
        constexpr double d =
            VOLUME_REF_TET * 10.0 / 189.0;
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

template<typename T>
void write_matrix(const T& matrix, const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (file.is_open()) {
        const uint64_t rows = matrix.rows();
        const uint64_t cols = matrix.cols();
        file.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
        file.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
        file.write(
            reinterpret_cast<const char*>(matrix.data()),
            rows * cols * sizeof(typename T::Scalar)
        );
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::_check_if_it_can_run() {
    _check_if_mesh_is_defined();
    if (!_is_any_source_defined){
        log::error("No sound source was defined."
            " Call add_sound_source to do so.");
    }
    if (!_is_freq_range_defined){
        log::error("Simulation frequency range was not defined."
            " Call set_frequency_range to do so.");
    }
    for (auto& evpg : _existing_evpg) {
        if (!_evpg_to_volprop.contains(evpg)) {
            log::error("Volume tag {} was not assigned."
            " Call add_volume_material to do so.", evpg);
        }
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::_define_freq_vector() {
    // todo: decide number here
    // todo: make it not linear
    _freq_steps = linspace(_freq_min, _freq_max, 1000);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::_organize_physical_group_data()
{
    // generate structures accessed through IPG
    _ispgv_to_volvel =
        fz::SafePtr<_FuncRealToCmplx>(_espg_to_volvel.size());
    for (const auto& [espg, volvel] : _espg_to_volvel) {
        const _ipg_t ispgv = _espg_to_ispg.at(espg);
        _ispgv_to_volvel[ispgv] = volvel;
    }
    _ispgp_to_pressure =
        fz::SafePtr<_FuncRealToCmplx>(_espg_to_pressure.size());
    for (const auto& [espg, pressure] : _espg_to_pressure) {
        const _ipg_t ispgp = _espg_to_ispg.at(espg);
        _ispgp_to_pressure[ispgp] = pressure;
    }
    _ispgi_to_impedance = fz::SafePtr<_FuncRealToCmplx>(_ispgi_count());
    for (const auto& [espg, impedance] : _espg_to_impedance) {
        const _ipg_t ispgi = _espg_to_ispg.at(espg);
        _ispgi_to_impedance[ispgi] = impedance;
    }
    _ivpg_to_volprop = fz::SafePtr<_VolProp>(_ivpg_count());
    for (const auto& [evpg, volprop] : _evpg_to_volprop) {
        const _ipg_t ivpg = _evpg_to_ivpg.at(evpg);
        _ivpg_to_volprop[ivpg] = volprop;
    }
    
    // generate the contiguous vector with the ispgi of each surface element
    size_t isei_count = 0;
    for (_idx_t sei=0; sei!=_sei_count(); ++sei) {
        if (_espg_to_impedance.contains(_sei_to_espg[sei])) {
            ++isei_count;
        }
    }
    _isei_to_sei = fz::SafePtr<_idx_t>(isei_count);
    isei_count = 0;
    for (_idx_t sei=0; sei!=_sei_count(); ++sei) {
        if (_espg_to_impedance.contains(_sei_to_espg[sei])) {
            const _idx_t new_isei = isei_count;
            _isei_to_sei[new_isei] = sei;
            ++isei_count;
        }
    }
    _isei_to_ispgi = fz::SafePtr<_ipg_t>(isei_count);
    for (_idx_t isei=0; isei!=_isei_count(); ++isei) {
        _idx_t sei = _isei_to_sei[isei];
        _isei_to_ispgi[isei] = _espg_to_ispg.at(_sei_to_espg[sei]);
    }
    
    // generate the contiguous vector with the ivpg of each volume element
    _vei_to_ivpg = fz::SafePtr<_ipg_t>(_vei_count());
    for (_idx_t vei=0; vei!=_vei_count(); ++vei) {
        _vei_to_ivpg[vei] = _evpg_to_ivpg.at(_vei_to_evpg[vei]);
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::_analyze_sparsity()
{
    std::unordered_set<std::pair<_idx_t,_idx_t>> existing_pairs;
    constexpr std::array<
        std::array<size_t,2>, COMB_REP_SIZE<NODES_IN_VOL_ELEM<O>,2>
    > COMBS_VOL = COMBINATION_REP<NODES_IN_VOL_ELEM<O>>;

    for (_idx_t vei=0; vei!=_vei_count(); ++vei) {
        for (const auto& c : COMBS_VOL) {
            existing_pairs.insert( 
                make_ascending_pair(
                    _vei_to_ni[vei][c[0]], _vei_to_ni[vei][c[1]]
                )
            );
        }
    }
    _nnz_rowcol_idx_pairs =
        fz::SafePtr<std::pair<_idx_t,_idx_t>>(existing_pairs.size());
    
    std::copy(
        existing_pairs.begin(),
        existing_pairs.end(),
        _nnz_rowcol_idx_pairs.begin()
    );
    std::sort(
        _nnz_rowcol_idx_pairs.begin(),
        _nnz_rowcol_idx_pairs.end(),
        compare_pair<_idx_t>
    );
}

template<ElementOrder O>
std::array<Eigen::Vector2d, NODES_IN_SFC_ELEM<O>> project_triangle_to_2d(
    const std::array<Eigen::Vector3d, NODES_IN_SFC_ELEM<O>>& vertices_3d
) {
    std::array<Eigen::Vector3d, NODES_IN_SFC_ELEM<O>> point_minus_origin;
    for (_idx_t ni = 0; ni!=NODES_IN_SFC_ELEM<O>; ++ni) {
        point_minus_origin[ni] = vertices_3d[ni] - vertices_3d[0];
    }
    const Eigen::Vector3d& u = point_minus_origin[1];
    const Eigen::Vector3d& v = point_minus_origin[2];
    const Eigen::Vector3d n = u.cross(v);
    const Eigen::Vector3d x = u / u.norm();
    Eigen::Vector3d y = n.cross(u);
    y /= y.norm();

    std::array<Eigen::Vector2d, NODES_IN_SFC_ELEM<O>> vertices_2d;
    for (_idx_t ni = 0; ni!=NODES_IN_SFC_ELEM<O>; ++ni) {
        vertices_2d[ni] = Eigen::Vector2d(
            point_minus_origin[ni].dot(x), point_minus_origin[ni].dot(y)
        );
    }
    return vertices_2d;
}

template<ElementOrder O>
void SimulationAcFemFreqD3<O>::_assemble_freq_independent_parts()
{   
    constexpr std::array<
        std::array<size_t,2>, COMB_REP_SIZE<NODES_IN_SFC_ELEM<O>,2>
    > COMBS_SFC = COMBINATION_REP<NODES_IN_SFC_ELEM<O>>;
    constexpr std::array<
        std::array<size_t,2>, COMB_REP_SIZE<NODES_IN_VOL_ELEM<O>,2>
    > COMBS_VOL = COMBINATION_REP<NODES_IN_VOL_ELEM<O>>;

    // count the nnz size for each ispgi
    fz::SafePtr<std::unordered_map<
        std::pair<_idx_t,_idx_t>, _idx_t
    >> ispgi_to_map_to_fipi(_ispgi_count());
    for (_idx_t isei=0; isei!=_isei_count(); ++isei)
    {
        const _ipg_t ispgi = _isei_to_ispgi[isei];
        const _idx_t sei = _isei_to_sei[isei];
        for (const auto& c : COMBS_SFC) {
            const std::pair<_idx_t,_idx_t> pair = make_ascending_pair(
                _sei_to_ni[sei][c[0]],
                _sei_to_ni[sei][c[1]]
            );
            if (!ispgi_to_map_to_fipi[ispgi].contains(pair)) {
                const _idx_t fipi = ispgi_to_map_to_fipi[ispgi].size();
                ispgi_to_map_to_fipi[ispgi].insert({pair, fipi});
            }
        }
    }

    // allocate memory in the safe pointers for the surface impedance elements
    _ispgi_to_damp_fi_part = fz::SafePtr<fz::SafePtr<_cmplx_t>>(_ispgi_count());
    _ispgi_to_ptr_in_a = fz::SafePtr<fz::SafePtr<_cmplx_t*>>(_ispgi_count());
    for (_ipg_t ispgi=0; ispgi!=_ispgi_count(); ++ispgi)
    {
        const size_t size = ispgi_to_map_to_fipi[ispgi].size();
        _ispgi_to_damp_fi_part[ispgi] = fz::SafePtr<_cmplx_t>(size);
        _ispgi_to_damp_fi_part[ispgi].fill(_cmplx_t(0.0, 0.0));
        _ispgi_to_ptr_in_a[ispgi] = fz::SafePtr<_cmplx_t*>(size);
    }

    // count the nnz size for each ivpg
    fz::SafePtr<std::unordered_map<
        std::pair<_idx_t,_idx_t>, _idx_t
    >> ivpg_to_map_to_fipi(_ivpg_count());
    for (_idx_t vei=0; vei!=_vei_count(); ++vei)
    {
        const _ipg_t ivpg = _vei_to_ivpg[vei];
        for (const auto& c : COMBS_VOL) {
            const std::pair<_idx_t,_idx_t> pair = make_ascending_pair(
                _vei_to_ni[vei][c[0]],
                _vei_to_ni[vei][c[1]]
            );
            if (!ivpg_to_map_to_fipi[ivpg].contains(pair)) {
                const _idx_t fipi = ivpg_to_map_to_fipi[ivpg].size();
                ivpg_to_map_to_fipi[ivpg].insert({pair, fipi});
            }
        }
    }

    // allocate memory in the safe pointers for the volume elements
    _ivpg_to_stif_fi_part = fz::SafePtr<fz::SafePtr<double>>(_ivpg_count());
    _ivpg_to_mass_fi_part = fz::SafePtr<fz::SafePtr<double>>(_ivpg_count());
    _ivpg_to_ptr_in_a = fz::SafePtr<fz::SafePtr<_cmplx_t*>>(_ivpg_count());
    for (_ipg_t ivpg=0; ivpg!=_ivpg_count(); ++ivpg)
    {
        const size_t size = ivpg_to_map_to_fipi[ivpg].size();
        _ivpg_to_stif_fi_part[ivpg] = fz::SafePtr<double>(size);
        _ivpg_to_stif_fi_part[ivpg].fill(0.0);
        _ivpg_to_mass_fi_part[ivpg] = fz::SafePtr<double>(size);
        _ivpg_to_mass_fi_part[ivpg].fill(0.0);
        _ivpg_to_ptr_in_a[ivpg] = fz::SafePtr<_cmplx_t*>(size);
    }
    
    // allocate the A matrix buffer
    _a_vals = fz::SafePtr<_cmplx_t>(_nnz_rowcol_idx_pairs.size());

    // loop over the surface impedance elements
    std::array<_idx_t, COMBS_SFC.size()> fipi_sfc;
    for (_idx_t isei=0; isei!=_isei_count(); ++isei)
    {
        const _ipg_t ispgi = _isei_to_ispgi[isei];
        const _ipg_t sei = _isei_to_sei[isei];
        
        // Create _ispgi_to_ptr_in_a
        for (_idx_t nci=0; nci!=COMBS_SFC.size(); ++nci)
        {
            const std::pair<_idx_t,_idx_t> pair = make_ascending_pair(
                _sei_to_ni[sei][COMBS_SFC[nci][0]],
                _sei_to_ni[sei][COMBS_SFC[nci][1]]
            );
            fipi_sfc[nci] = ispgi_to_map_to_fipi[ispgi].at(pair);
            
            const std::pair<_idx_t,_idx_t>* const pair_ptr = std::lower_bound(
                _nnz_rowcol_idx_pairs.begin(),
                _nnz_rowcol_idx_pairs.end(),
                pair,
                compare_pair<_idx_t>
            );
            const ptrdiff_t ptrdiff = pair_ptr - _nnz_rowcol_idx_pairs.begin();
            _ispgi_to_ptr_in_a[ispgi][fipi_sfc[nci]] =
                _a_vals.begin() + ptrdiff;
        }

        // coordinates matrix
        std::array<Eigen::Vector3d, NODES_IN_SFC_ELEM<O>> triangle_3d;
        for (size_t ni=0; ni!=NODES_IN_SFC_ELEM<O>; ++ni) {
            const _idx_t node_idx = _sei_to_ni[sei][ni];
            triangle_3d[ni] = Eigen::Vector3d(
                _node_coords[node_idx][0],
                _node_coords[node_idx][1],
                _node_coords[node_idx][2]
            );
        }
        std::array<Eigen::Vector2d, NODES_IN_SFC_ELEM<O>> triangle_2d =
            project_triangle_to_2d<O>(triangle_3d);
        Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,2> coords_matrix;
        for (size_t ni=0; ni!=NODES_IN_SFC_ELEM<O>; ++ni) {
            coords_matrix(ni,0) = triangle_2d[ni](0);
            coords_matrix(ni,1) = triangle_2d[ni](1);
        }

        // damping matrix
        constexpr std::array<std::array<double,2>,NGP_DAMP<O>>
            GAUSS_POINTS_DAMP = GAUSS_POINTS_SFC<NGP_DAMP<O>>;
        for (_idx_t gpi=0; gpi!=NGP_DAMP<O>; ++gpi)
        {
            Eigen::Matrix<double,2,NODES_IN_SFC_ELEM<O>> nabla_n =
                shape_func_sfc_gradient<O>(
                    GAUSS_POINTS_DAMP[gpi][0],
                    GAUSS_POINTS_DAMP[gpi][1]
                ); // todo: try putting constexpr here

            const Eigen::Matrix<double,2,2> jac_matrix =
                nabla_n * coords_matrix;
            
            const double det_jac = jac_matrix.determinant();
            
            // std::array<std::array<double,DIM>,3> triangle_coords;
            // for (size_t ni=0; ni!=3; ++ni) {
            //     const _idx_t node_idx = _sei_to_ni[sei][ni];
            //     triangle_coords[ni] = std::array<double,DIM>({
            //         _node_coords[node_idx][0],
            //         _node_coords[node_idx][1],
            //         _node_coords[node_idx][2]
            //     });
            // }
            // const double det_jac = 2*get_triangle_area(triangle_coords);

            const Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,1> n =
                shape_func_sfc<O>(
                    GAUSS_POINTS_DAMP[gpi][0],
                    GAUSS_POINTS_DAMP[gpi][1]
                );
            
            const
            Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,NODES_IN_SFC_ELEM<O>>
                nnt = n * n.transpose();

            // todo: multiply detj and w without creating another eigen matrix
            const
            Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,NODES_IN_SFC_ELEM<O>>
                nnt_detj_w =
                    nnt * det_jac * GAUSS_WEIGHTS_SFC<NGP_DAMP<O>>[gpi];
            
            for (_idx_t nci=0; nci!=COMBS_SFC.size(); ++nci) {
                _ispgi_to_damp_fi_part[ispgi][fipi_sfc[nci]] += nnt_detj_w(
                    COMBS_SFC[nci][0], COMBS_SFC[nci][1]
                );
            }
        }
    }

    // loop over the volume elements
    std::array<_idx_t, COMBS_VOL.size()> fipi_vol;
    for (_idx_t vei=0; vei!=_vei_count(); ++vei)
    {
        const _ipg_t ivpg = _vei_to_ivpg[vei];
        
        // Create _ivpg_to_ptr_in_a
        for (_idx_t nci=0; nci!=COMBS_VOL.size(); ++nci)
        {
            const std::pair<_idx_t,_idx_t> pair = make_ascending_pair(
                _vei_to_ni[vei][COMBS_VOL[nci][0]],
                _vei_to_ni[vei][COMBS_VOL[nci][1]]
            );
            fipi_vol[nci] = ivpg_to_map_to_fipi[ivpg].at(pair);
            
            const std::pair<_idx_t,_idx_t>* const pair_ptr = std::lower_bound(
                _nnz_rowcol_idx_pairs.begin(),
                _nnz_rowcol_idx_pairs.end(),
                pair,
                compare_pair<_idx_t>
            );
            const ptrdiff_t ptrdiff = pair_ptr - _nnz_rowcol_idx_pairs.begin();
            _ivpg_to_ptr_in_a[ivpg][fipi_vol[nci]] = _a_vals.begin() + ptrdiff;
        }
        
        // coordinates matrix
        Eigen::Matrix<double,NODES_IN_VOL_ELEM<O>,DIM> coords_matrix;
        for (size_t ni=0; ni!=NODES_IN_VOL_ELEM<O>; ++ni) {
            const _idx_t node_idx = _vei_to_ni[vei][ni];
            for (size_t di=0; di!=DIM; ++di) {
                coords_matrix(ni,di) = _node_coords[node_idx][di];
            }
        }

        // stiffness matrix
        constexpr std::array<std::array<double,DIM>,NGP_STIF<O>>
            GAUSS_POINTS_STIF = GAUSS_POINTS_VOL<NGP_STIF<O>>;
        for (_idx_t gpi=0; gpi!=NGP_STIF<O>; ++gpi)
        {
            Eigen::Matrix<double,DIM,NODES_IN_VOL_ELEM<O>> nabla_n =
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

            for (_idx_t nci=0; nci!=COMBS_VOL.size(); ++nci) {
                _ivpg_to_stif_fi_part[ivpg][fipi_vol[nci]] += btb_detj_w(
                    COMBS_VOL[nci][0], COMBS_VOL[nci][1]
                );
            }
        }

        // mass matrix
        constexpr std::array<std::array<double,DIM>,NGP_MASS<O>>
            GAUSS_POINTS_MASS = GAUSS_POINTS_VOL<NGP_MASS<O>>;
        for (_idx_t gpi=0; gpi!=NGP_MASS<O>; ++gpi)
        {
            Eigen::Matrix<double,DIM,NODES_IN_VOL_ELEM<O>> nabla_n =
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
            
            for (_idx_t nci=0; nci!=COMBS_VOL.size(); ++nci) {
                _ivpg_to_mass_fi_part[ivpg][fipi_vol[nci]] += nnt_detj_w(
                    COMBS_VOL[nci][0], COMBS_VOL[nci][1]
                );
            }
        }
    }

    ispgi_to_map_to_fipi.free();
    ivpg_to_map_to_fipi.free();
}

void solve_using_onemkl(
    const fz::SafePtr<_cmplx_t>& a_vals,
    const fz::SafePtr<std::pair<_idx_t,_idx_t>>& nnz_rowcol_idx_pairs,
    const fz::SafePtr<_cmplx_t>& b,
    fz::SafePtr<_cmplx_t>& x
) {
    auto print_dss_error = [](const _INTEGER_t* const error_id) {
        fprintf(stderr, "MLK code %lli\n", *error_id);
        exit(1);
        return;
    };

    // problem dimensions
    const MKL_INT node_count = b.size();
    const MKL_INT nnz_count = a_vals.size();

    // define sparsity patterns
    fz::SafePtr<MKL_INT> a_col_idx(nnz_count);
    fz::SafePtr<MKL_INT> a_row_ptr(node_count + 1);
    MKL_INT* it_a_col_idx = a_col_idx.begin();
    MKL_INT* it_a_row_ptr = a_row_ptr.begin();
    const std::pair<_idx_t,_idx_t>* it_nnz_rowcol_idx_pairs =
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
    
    // system properties
    const MKL_INT symmetry_type = MKL_DSS_SYMMETRIC_COMPLEX;
    const MKL_INT positive_definiteness = MKL_DSS_INDEFINITE;

    // options
    const MKL_INT options = 
        MKL_DSS_MSG_LVL_WARNING +
        MKL_DSS_TERM_LVL_ERROR +
        MKL_DSS_ZERO_BASED_INDEXING;

    // allocate memory for the solver handle and error id.
    _MKL_DSS_HANDLE_t dss_handle;
    _INTEGER_t error_id;
    
    // initialize the solver
    error_id = dss_create(dss_handle, options);
    if (error_id != MKL_DSS_SUCCESS) { print_dss_error(&error_id); }

    // define the non-zero structure of the matrix
    error_id = dss_define_structure(
        dss_handle, symmetry_type, a_row_ptr.data(),
        node_count, node_count, a_col_idx.data(), nnz_count
    );
    if (error_id != MKL_DSS_SUCCESS) { print_dss_error(&error_id); }
    
    // reorder the matrix
    error_id = dss_reorder(dss_handle, options, 0);
    if (error_id != MKL_DSS_SUCCESS) { print_dss_error(&error_id); }
    
    // factor the matrix
    error_id = dss_factor_complex(
        dss_handle, positive_definiteness,
        reinterpret_cast<const double*>(a_vals.data())
    );
    if (error_id != MKL_DSS_SUCCESS) { print_dss_error(&error_id); }
    
    // solution vector
    const MKL_INT num_of_b = 1;
    error_id = dss_solve_complex(
        dss_handle, options,
        reinterpret_cast<const double*>(b.data()),
        num_of_b, reinterpret_cast<double*>(x.data())
    );
    if (error_id != MKL_DSS_SUCCESS) { print_dss_error(&error_id); }
    a_row_ptr.free();
    a_col_idx.free();
}

void solve_using_eigen(
    const fz::SafePtr<_cmplx_t>& a_vals,
    const fz::SafePtr<std::pair<_idx_t,_idx_t>>& nnz_rowcol_idx_pairs,
    const fz::SafePtr<_cmplx_t>& b,
    fz::SafePtr<_cmplx_t>& x_out
) {
    const size_t node_count = b.size();
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
    Eigen::Matrix<_cmplx_t, Eigen::Dynamic, 1> b_eig(node_count);
    std::copy(b.begin(), b.end(), b_eig.data());

    Eigen::SparseLU<Eigen::SparseMatrix<_cmplx_t>> solver;
    solver.analyzePattern(a);
    solver.factorize(a);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Factorization failed. Matrix may be singular.\n";
    }
    const Eigen::VectorXcd x = solver.solve(b_eig);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solving failed.\n";
    }
    for (size_t j=0; j!=node_count; ++j) {
        x_out[j] = x(j);
    }
}

template<ElementOrder O>
void SimulationAcFemFreqD3<O>::_solve()
{
    _result = ResultAcFemFreqD3(_ni_count(), _freq_count());
    fz::SafePtr<_cmplx_t> b(_ni_count());

    // start timer
    auto start_time = std::chrono::system_clock::now();
    auto time_t_start = std::chrono::system_clock::to_time_t(start_time);
    std::cout << "Solver started at: "
        << std::put_time(std::localtime(&time_t_start), "%H:%M:%S") << "\n";
    std::cout << "Solution progress: 0%\n";

    for (size_t fi=0; fi!=_freq_count(); ++fi)
    {
        _a_vals.fill(_cmplx_t(0.0, 0.0));
        const double freq = _freq_steps[fi];
        const double omega = 2*std::numbers::pi*freq;
        const double omega_squared = std::pow(omega,2);

        // add damping matrix to a
        for (_idx_t ispgi=0; ispgi!=_ispgi_count(); ++ispgi) {
            const _cmplx_t impedance_value = _ispgi_to_impedance[ispgi](freq);
            const _cmplx_t damp_fd_part =
                std::complex<double>(0.0, omega) / impedance_value;
            
            for (_idx_t fipi=0; fipi!=_ispgi_to_ptr_in_a[ispgi].size(); ++fipi)
            {
                *_ispgi_to_ptr_in_a[ispgi][fipi] +=
                    _ispgi_to_damp_fi_part[ispgi][fipi] * damp_fd_part;
            }
        }

        // add stiffness and mass matrix to a
        for (_idx_t ivpg=0; ivpg!=_ivpg_count(); ++ivpg) {
            const _cmplx_t density_value =
                (_ivpg_to_volprop[ivpg].density)(freq);
            const _cmplx_t soundspeed_value =
                (_ivpg_to_volprop[ivpg].soundspeed)(freq);

            const _cmplx_t stif_fd_part = 1.0 / density_value;
            const _cmplx_t mass_fd_part =
                - omega_squared / (density_value*std::pow(soundspeed_value,2));
            
            for (_idx_t fipi=0; fipi!=_ivpg_to_ptr_in_a[ivpg].size(); ++fipi) {
                *_ivpg_to_ptr_in_a[ivpg][fipi] +=
                    _ivpg_to_stif_fi_part[ivpg][fipi] * stif_fd_part +
                    _ivpg_to_mass_fi_part[ivpg][fipi] * mass_fd_part;
            }
        }

        // b vector
        b.fill(0);
        for (size_t j=0; j!=_point_volvel.size(); ++j)
        {
            const _cmplx_t volvel =
                std::get<_FuncRealToCmplx>(_point_volvel[j])(freq);

            b[std::get<_idx_t>(_point_volvel[j])] =
                _cmplx_t(0.0,-1.0)*omega*volvel;
        }

        // solve
        fz::SafePtr<_cmplx_t> x(_ni_count());
        #if NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
            solve_using_eigen(_a_vals, _nnz_rowcol_idx_pairs, b, x);
        #elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
            solve_using_onemkl(_a_vals, _nnz_rowcol_idx_pairs, b, x);
        #else
            static_assert(false, "Invalid NUMAV_SYSTEM_SOLVER.");
        #endif
        for (_idx_t n=0; n!=_ni_count(); ++n) {
            _result._data(n,fi) = x[n];
        }
        for (_idx_t n=0; n!=_ni_count(); ++n) {
            assert(_result._data(n,fi) == x[n]);
        }
        x.free();

        // print the progress
        std::cout << "\033[A"; // Moves the cursor up one line
        std::cout << "\033[2K"; // Clears the entire line
        double progress_percentage =
            100.0*static_cast<double>(fi+1)/static_cast<double>(_freq_count());
        std::cout << "Solution progress: " << std::fixed << 
            std::setprecision(2) << progress_percentage << "%";

        if (fi > 99) { // print the time left
            auto current_time = std::chrono::system_clock::now();
            auto time_taken_until_now =
                std::chrono::duration_cast<std::chrono::seconds>(
                    current_time - start_time
                );
            auto time_left = 
                time_taken_until_now.count()/(fi+1)*(_freq_count()-(fi+1));
            
            auto hours = time_left / 3600;
            auto minutes = (time_left % 3600) / 60;
            std::cout << " - estimated time left: " <<
                hours << " h and " << minutes << " min\n";
        }
        else {
            std::cout << "\n";
        }
    }
    b.free();
    write_matrix(_result._data, "pressure.bin");

    // print finish
    auto end_time = std::chrono::system_clock::now();
    auto time_t_end = std::chrono::system_clock::to_time_t(end_time);
    std::cout << "Solver ended at: " << std::put_time(std::localtime(&time_t_end), "%H:%M:%S") << "\n";
}

template <ElementOrder O>
ResultAcFemFreqD3 SimulationAcFemFreqD3<O>::run()
{
    _check_if_it_can_run();
    log::print_opening();
    log::print_opening_ac_fem_freq_d3();
    _define_freq_vector();
    _organize_physical_group_data();
    _analyze_sparsity();
    _assemble_freq_independent_parts();
    _solve();
    return ResultAcFemFreqD3();
}

} // namespace numav
