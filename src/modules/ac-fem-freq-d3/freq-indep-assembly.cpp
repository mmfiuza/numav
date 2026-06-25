// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/impl.hpp"

#include "common/exception.hpp"
#include "common/maths.hpp"
#include "common/hash-functions.hpp"
#include "common/utils.hpp"
#include "modules/ac-fem-freq-d3/shape-functions.hpp"
#include "modules/ac-fem-freq-d3/analytic-integration.hpp"
#include "modules/ac-fem-freq-d3/gauss-quadrature.hpp"

#include <set>

namespace numav {

template<typename T>
bool compare_pair(const std::pair<T,T> a, const std::pair<T,T> b) {
    if constexpr
    (GLOBAL_MATRIX_STORAGE_ORDER == MatrixStorageOrder::COL_MAJOR)
    {
        if (a.second != b.second) {
            return a.second < b.second;
        }
        else {
            return a.first < b.first;
        }
    }
    else if constexpr
    (GLOBAL_MATRIX_STORAGE_ORDER == MatrixStorageOrder::ROW_MAJOR)
    {
        if (a.first != b.first) {
            return a.first < b.first;
        }
        else {
            return a.second < b.second;
        }
    }
}

template<typename T>
std::pair<T,T> make_ordered_rowcol_pair(const T a, const T b)
{
    #if NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
        return std::make_pair(a,b);
    #endif

    if constexpr
    (GLOBAL_MATRIX_TRIANGULAR_TYPE == TriangularMatrixType::UPPER)
    {
        return a<b ? std::make_pair(a,b) : std::make_pair(b,a);
    }
    else if constexpr
    (GLOBAL_MATRIX_TRIANGULAR_TYPE == TriangularMatrixType::LOWER)
    {
        return a<b ? std::make_pair(b,a) : std::make_pair(a,b);
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_allocate_a()
{
    #if NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
        constexpr std::array<
            std::array<size_t, 2UL>,
            PERMUTATION_REP_SIZE<ENI_COUNT_VOL<O>, 2UL>
        > ENI_PAIRS = PERMUTATION_REP<ENI_COUNT_VOL<O>>;
    #else
        constexpr std::array<
            std::array<size_t, 2UL>,
            COMBINATION_REP_SIZE<ENI_COUNT_VOL<O>, 2UL>
        > ENI_PAIRS = COMBINATION_REP<ENI_COUNT_VOL<O>>;
    #endif
    std::unordered_set<std::pair<size_t,size_t>> existing_pairs;
    for (size_t vei = 0UL; vei != _vei_count; ++vei) {
        for (const auto& eni : ENI_PAIRS) {
            #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
                if (eni[0UL] == eni[1UL]) { continue; }
            #endif
            existing_pairs.insert(
                make_ordered_rowcol_pair(
                    _vei_to_ni[vei][eni[0UL]], _vei_to_ni[vei][eni[1UL]]
                )
            );
        }
    }
    _ni_connections = fz::SafePtr<std::pair<size_t,size_t>>(
        existing_pairs.size()
    );
    std::copy(
        existing_pairs.begin(),
        existing_pairs.end(),
        _ni_connections.begin()
    );
    std::sort(
        _ni_connections.begin(),
        _ni_connections.end(),
        compare_pair<size_t>
    );
    _a_vals = fz::SafePtr<Cmplx>(_ni_connections.size());
    #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
        _a_diag = fz::SafePtr<Cmplx>(_ni_count);
    #endif
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_allocate_b()
{
    std::unordered_set<size_t> existing_source_nodes;
    for (size_t vsei = 0UL; vsei != _vsei_count; ++vsei) {
        const size_t sei = _vsei_to_sei[vsei];
        for (size_t eni = 0UL; eni != ENI_COUNT_SFC<O>; ++eni) {
            const size_t ni = _sei_to_ni[sei][eni];
            existing_source_nodes.insert(ni);
        }
    }
    for (size_t vpi = 0UL; vpi != _vpi_count; ++vpi) {
        existing_source_nodes.insert(_vpi_to_ni[vpi]);
    }
    for (size_t psei = 0UL; psei != _psei_count; ++psei) {
        const size_t sei = _psei_to_sei[psei];
        for (size_t eni = 0UL; eni != ENI_COUNT_SFC<O>; ++eni) {
            const size_t ni = _sei_to_ni[sei][eni];
            existing_source_nodes.insert(ni);
        }
    }
    for (size_t ppi = 0UL; ppi != _ppi_count; ++ppi) {
        existing_source_nodes.insert(_ppi_to_ni[ppi]);
    }
    _b_row_idx = fz::SafePtr<size_t>(existing_source_nodes.size());
    std::copy(
        existing_source_nodes.begin(),
        existing_source_nodes.end(),
        _b_row_idx.begin()
    );
    std::sort(_b_row_idx.begin(), _b_row_idx.end());
    _b_vals = fz::SafePtr<Cmplx>(_b_row_idx.size());
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_allocate_x()
{
    _x = fz::SafePtr<Cmplx>(_ni_count);
}

template<ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_assemble_fi_part_for_vol_elements()
{
    #if NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
        constexpr std::array<
            std::array<size_t, 2UL>,
            PERMUTATION_REP_SIZE<ENI_COUNT_VOL<O>, 2UL>
        > ENI_PAIRS = PERMUTATION_REP<ENI_COUNT_VOL<O>>;
    #else
        constexpr std::array<
            std::array<size_t, 2UL>,
            COMBINATION_REP_SIZE<ENI_COUNT_VOL<O>, 2UL>
        > ENI_PAIRS = COMBINATION_REP<ENI_COUNT_VOL<O>>;
    #endif

    // count the fipi for each ivpg
    fz::SafePtr<std::unordered_map<
        std::pair<size_t,size_t>, size_t
    >> ivpg_to_map_to_fipi(_ivpg_count);
    for (size_t vei = 0UL; vei != _vei_count; ++vei)
    {
        const size_t ivpg = _vei_to_ivpg[vei];
        for (const auto& eni : ENI_PAIRS) {
            const std::pair<size_t,size_t> pair = make_ordered_rowcol_pair(
                _vei_to_ni[vei][eni[0UL]], _vei_to_ni[vei][eni[1UL]]
            );
            if (!(ivpg_to_map_to_fipi[ivpg].count(pair) > 0)) {
                const size_t fipi = ivpg_to_map_to_fipi[ivpg].size();
                ivpg_to_map_to_fipi[ivpg].insert({pair, fipi});
            }
        }
    }

    // allocate memory in the safe pointers
    _ivpg_to_stif_fi_part = fz::SafePtr<fz::SafePtr<Float>>(_ivpg_count);
    _ivpg_to_mass_fi_part = fz::SafePtr<fz::SafePtr<Float>>(_ivpg_count);
    _ivpg_to_ptr_in_a = fz::SafePtr<fz::SafePtr<Cmplx*>>(_ivpg_count);
    for (size_t ivpg = 0UL; ivpg != _ivpg_count; ++ivpg)
    {
        const size_t size = ivpg_to_map_to_fipi[ivpg].size();
        _ivpg_to_stif_fi_part[ivpg] = fz::SafePtr<Float>(size);
        _ivpg_to_stif_fi_part[ivpg].fill(0_F);
        _ivpg_to_mass_fi_part[ivpg] = fz::SafePtr<Float>(size);
        _ivpg_to_mass_fi_part[ivpg].fill(0_F);
        _ivpg_to_ptr_in_a[ivpg] = fz::SafePtr<Cmplx*>(size);
    }

    // assemble elemental stiffness and mass matrices
    std::array<size_t, ENI_PAIRS.size()> fipi_vol;
    for (size_t vei = 0UL; vei != _vei_count; ++vei)
    {
        const size_t ivpg = _vei_to_ivpg[vei];
        
        // Create _ivpg_to_ptr_in_a
        for (size_t nci = 0UL; nci != ENI_PAIRS.size(); ++nci)
        {
            const std::pair<size_t,size_t> pair = make_ordered_rowcol_pair(
                _vei_to_ni[vei][ENI_PAIRS[nci][0UL]],
                _vei_to_ni[vei][ENI_PAIRS[nci][1UL]]
            );
            fipi_vol[nci] = ivpg_to_map_to_fipi[ivpg].at(pair);

            #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
                if(pair.first == pair.second) {
                    const size_t ni = pair.first;
                    _ivpg_to_ptr_in_a[ivpg][fipi_vol[nci]] =
                        _a_diag.data() + ni;
                    continue;
                }
            #endif

            const std::pair<size_t,size_t>* const pair_ptr = std::lower_bound(
                _ni_connections.begin(),
                _ni_connections.end(),
                pair,
                compare_pair<size_t>
            );
            const size_t ptrdiff = pair_ptr - _ni_connections.begin();
            _ivpg_to_ptr_in_a[ivpg][fipi_vol[nci]] = _a_vals.begin() + ptrdiff;
        }
        
        // coordinates matrix
        Eigen::Matrix<Float,DIM,ENI_COUNT_VOL<O>> coords_matrix;
        for (size_t eni = 0UL; eni != ENI_COUNT_VOL<O>; ++eni) {
            const size_t ni = _vei_to_ni[vei][eni];
            for (size_t di = 0UL; di != DIM; ++di) {
                coords_matrix(di,eni) = _ni_to_coords[ni][di];
            }
        }
        #if \
        NUMAV_STIF_INTEGRATION_METHOD == NUMAV_ANALYTIC || \
        NUMAV_MASS_INTEGRATION_METHOD == NUMAV_ANALYTIC
            const Float tet_volume = get_tetrahedron_volume(
                *reinterpret_cast<std::array<std::array<Float,3UL>,4UL>*>(
                    coords_matrix.data()
                )
            );
        #endif

        // stiffness matrix
        #if NUMAV_STIF_INTEGRATION_METHOD == NUMAV_GAUSS_QUADRATURE
            constexpr std::array<std::array<Float,DIM>,NGP_STIF<O>>
                GAUSS_POINTS_STIF = GAUSS_POINTS_VOL<NGP_STIF<O>>;
            for (size_t gpi = 0UL; gpi != NGP_STIF<O>; ++gpi)
            {
                const Eigen::Matrix<Float,ENI_COUNT_VOL<O>,DIM> nabla_n =
                    shape_func_vol_gradient<O>(
                        GAUSS_POINTS_STIF[gpi][0UL],
                        GAUSS_POINTS_STIF[gpi][1UL],
                        GAUSS_POINTS_STIF[gpi][2UL]
                    ); // todo: try putting constexpr here

                const Eigen::Matrix<Float,DIM,DIM> jac_matrix = 
                    coords_matrix * nabla_n;

                const Eigen::Matrix<Float,DIM,DIM> inv_jac =
                    jac_matrix.inverse();
                
                const Eigen::Matrix<Float,ENI_COUNT_VOL<O>,DIM> b_matrix =
                    nabla_n * inv_jac;
                
                const Eigen::Matrix<Float,ENI_COUNT_VOL<O>,ENI_COUNT_VOL<O>>
                    bbt = b_matrix * b_matrix.transpose();
                
                const Float det_jac = std::abs(jac_matrix.determinant());
                
                // todo: multiply detj and w without creating another matrix
                const Eigen::Matrix<Float,ENI_COUNT_VOL<O>,ENI_COUNT_VOL<O>>
                    bbt_detj_w = 
                        bbt * det_jac * GAUSS_WEIGHTS_VOL<NGP_STIF<O>>[gpi];

                for (size_t nci = 0UL; nci != ENI_PAIRS.size(); ++nci) {
                    _ivpg_to_stif_fi_part[ivpg][fipi_vol[nci]] +=
                        bbt_detj_w(ENI_PAIRS[nci][0UL], ENI_PAIRS[nci][1UL]);
                }
            }
        #elif NUMAV_STIF_INTEGRATION_METHOD == NUMAV_ANALYTIC
            const Eigen::Matrix<Float, ENI_COUNT_VOL<O>, ENI_COUNT_VOL<O>>
                stif_matrix_const_part = 
                    get_stif_matrix_const_part<O>(coords_matrix);

            const Eigen::Matrix<Float, ENI_COUNT_VOL<O>, ENI_COUNT_VOL<O>>
                stif_fi_part = stif_matrix_const_part / tet_volume;

            for (size_t nci = 0UL; nci != ENI_PAIRS.size(); ++nci) {
                _ivpg_to_stif_fi_part[ivpg][fipi_vol[nci]] +=
                    stif_fi_part(ENI_PAIRS[nci][0UL], ENI_PAIRS[nci][1UL]);
            }
        #endif

        // mass matrix
        #if NUMAV_MASS_INTEGRATION_METHOD == NUMAV_GAUSS_QUADRATURE
            constexpr std::array<std::array<Float,DIM>,NGP_MASS<O>>
                GAUSS_POINTS_MASS = GAUSS_POINTS_VOL<NGP_MASS<O>>;
            for (size_t gpi = 0UL; gpi != NGP_MASS<O>; ++gpi)
            {
                const Eigen::Matrix<Float,ENI_COUNT_VOL<O>,DIM> nabla_n =
                    shape_func_vol_gradient<O>(
                        GAUSS_POINTS_MASS[gpi][0UL],
                        GAUSS_POINTS_MASS[gpi][1UL],
                        GAUSS_POINTS_MASS[gpi][2UL]
                    ); // todo: try putting constexpr here

                const Eigen::Matrix<Float,DIM,DIM> jac_matrix = 
                    coords_matrix * nabla_n;
                    
                const Eigen::Matrix<Float,ENI_COUNT_VOL<O>,1UL> n =
                    shape_func_vol<O>(
                        GAUSS_POINTS_MASS[gpi][0UL],
                        GAUSS_POINTS_MASS[gpi][1UL],
                        GAUSS_POINTS_MASS[gpi][2UL]
                    );
                    
                const Eigen::Matrix<Float,ENI_COUNT_VOL<O>,ENI_COUNT_VOL<O>>
                    nnt = n * n.transpose();

                const Float det_jac = std::abs(jac_matrix.determinant());
                    
                // todo: multiply detj and w without creating another matrix
                const Eigen::Matrix<Float,ENI_COUNT_VOL<O>,ENI_COUNT_VOL<O>>
                    nnt_detj_w =
                        nnt * det_jac * GAUSS_WEIGHTS_VOL<NGP_MASS<O>>[gpi];
                
                for (size_t nci = 0UL; nci != ENI_PAIRS.size(); ++nci) {
                    _ivpg_to_mass_fi_part[ivpg][fipi_vol[nci]] +=
                        nnt_detj_w(ENI_PAIRS[nci][0UL], ENI_PAIRS[nci][1UL]);
                }
            }
        #elif NUMAV_MASS_INTEGRATION_METHOD == NUMAV_ANALYTIC
            const Eigen::Matrix<Float,ENI_COUNT_VOL<O>,ENI_COUNT_VOL<O>>
                mass_matrix_const_part = MASS_MATRIX_CONST_PART<O>;
            
            const Eigen::Matrix<Float,ENI_COUNT_VOL<O>,ENI_COUNT_VOL<O>>
                mass_fi_part = mass_matrix_const_part * tet_volume;
            
            for (size_t nci = 0UL; nci != ENI_PAIRS.size(); ++nci) {
                _ivpg_to_mass_fi_part[ivpg][fipi_vol[nci]] +=
                    mass_fi_part(ENI_PAIRS[nci][0UL], ENI_PAIRS[nci][1UL]);
            }
        #endif
    }
    ivpg_to_map_to_fipi.free();
}

template<ElementOrder O>
std::array<Eigen::Vector2d, ENI_COUNT_SFC<O>> project_triangle_to_2d(
    const std::array<Eigen::Vector3d, ENI_COUNT_SFC<O>> vertices_3d
) {
    std::array<Eigen::Vector3d, ENI_COUNT_SFC<O>> point_minus_origin;
    for (size_t eni = 0UL; eni != ENI_COUNT_SFC<O>; ++eni) {
        point_minus_origin[eni] = vertices_3d[eni] - vertices_3d[0UL];
    }
    const Eigen::Vector3d& u = point_minus_origin[1UL];
    const Eigen::Vector3d& v = point_minus_origin[2UL];
    const Eigen::Vector3d n = u.cross(v);
    const Eigen::Vector3d x = u / u.norm();
    Eigen::Vector3d y = n.cross(u);
    y /= y.norm();

    std::array<Eigen::Vector2d, ENI_COUNT_SFC<O>> vertices_2d;
    for (size_t eni = 0UL; eni != ENI_COUNT_SFC<O>; ++eni) {
        vertices_2d[eni] = Eigen::Vector2d(
            point_minus_origin[eni].dot(x), point_minus_origin[eni].dot(y)
        );
    }
    return vertices_2d;
}

template<ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_assemble_fi_part_for_sfc_impedance()
{
    #if NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
        constexpr std::array<
            std::array<size_t, 2UL>,
            PERMUTATION_REP_SIZE<ENI_COUNT_SFC<O>, 2UL>
        > ENI_PAIRS = PERMUTATION_REP<ENI_COUNT_SFC<O>>;
    #else
        constexpr std::array<
            std::array<size_t, 2UL>,
            COMBINATION_REP_SIZE<ENI_COUNT_SFC<O>, 2UL>
        > ENI_PAIRS = COMBINATION_REP<ENI_COUNT_SFC<O>>;
    #endif

    // count the fipi for each ispgi
    fz::SafePtr<std::unordered_map<
        std::pair<size_t,size_t>, size_t
    >> ispgi_to_map_to_fipi(_ispgi_count);
    for (size_t isei = 0UL; isei != _isei_count; ++isei)
    {
        const size_t ispgi = _isei_to_ispgi[isei];
        const size_t sei = _isei_to_sei[isei];
        for (const auto& eni : ENI_PAIRS) {
            const std::pair<size_t,size_t> pair = make_ordered_rowcol_pair(
                _sei_to_ni[sei][eni[0UL]], _sei_to_ni[sei][eni[1UL]]
            );
            if (!(ispgi_to_map_to_fipi[ispgi].count(pair) > 0)) {
                const size_t fipi = ispgi_to_map_to_fipi[ispgi].size();
                ispgi_to_map_to_fipi[ispgi].insert({pair, fipi});
            }
        }
    }

    // allocate memory in the safe pointers
    _ispgi_to_damp_fi_part = fz::SafePtr<fz::SafePtr<Float>>(_ispgi_count);
    _ispgi_to_ptr_in_a = fz::SafePtr<fz::SafePtr<Cmplx*>>(_ispgi_count);
    for (size_t ispgi = 0UL; ispgi != _ispgi_count; ++ispgi)
    {
        const size_t size = ispgi_to_map_to_fipi[ispgi].size();
        _ispgi_to_damp_fi_part[ispgi] = fz::SafePtr<Float>(size);
        _ispgi_to_damp_fi_part[ispgi].fill(0_F);
        _ispgi_to_ptr_in_a[ispgi] = fz::SafePtr<Cmplx*>(size);
    }

    // assemble the elemental damping matrices
    std::array<size_t, ENI_PAIRS.size()> fipi_damp;
    for (size_t isei = 0UL; isei != _isei_count; ++isei)
    {
        const size_t ispgi = _isei_to_ispgi[isei];
        const size_t sei = _isei_to_sei[isei];
        
        // Create _ispgi_to_ptr_in_a
        for (size_t nci = 0UL; nci != ENI_PAIRS.size(); ++nci)
        {
            const std::pair<size_t,size_t> pair = make_ordered_rowcol_pair(
                _sei_to_ni[sei][ENI_PAIRS[nci][0UL]],
                _sei_to_ni[sei][ENI_PAIRS[nci][1UL]]
            );
            fipi_damp[nci] = ispgi_to_map_to_fipi[ispgi].at(pair);
            
            #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
                if(pair.first == pair.second) {
                    const size_t ni = pair.first;
                    _ispgi_to_ptr_in_a[ispgi][fipi_damp[nci]] =
                        _a_diag.data() + ni;
                    continue;
                }
            #endif

            const std::pair<size_t,size_t>* const pair_ptr = std::lower_bound(
                _ni_connections.begin(),
                _ni_connections.end(),
                pair,
                compare_pair<size_t>
            );
            const size_t ptrdiff = pair_ptr - _ni_connections.begin();
            _ispgi_to_ptr_in_a[ispgi][fipi_damp[nci]] =
                _a_vals.begin() + ptrdiff;
        }

        #if NUMAV_DAMP_INTEGRATION_METHOD == NUMAV_GAUSS_QUADRATURE
            // coordinates matrix
            std::array<Eigen::Vector3d, ENI_COUNT_SFC<O>> triangle_3d;
            for (size_t eni = 0UL; eni != ENI_COUNT_SFC<O>; ++eni) {
                const size_t ni = _sei_to_ni[sei][eni];
                triangle_3d[eni] = Eigen::Vector3d(
                    _ni_to_coords[ni][0UL],
                    _ni_to_coords[ni][1UL],
                    _ni_to_coords[ni][2UL]
                );
            }
            std::array<Eigen::Vector2d, ENI_COUNT_SFC<O>> triangle_2d =
                project_triangle_to_2d<O>(triangle_3d);
            Eigen::Matrix<Float, 2UL, ENI_COUNT_SFC<O>> coords_matrix;
            for (size_t eni = 0UL; eni != ENI_COUNT_SFC<O>; ++eni) {
                coords_matrix(0UL,eni) = triangle_2d[eni](0UL);
                coords_matrix(1UL,eni) = triangle_2d[eni](1UL);
            }
        #elif NUMAV_DAMP_INTEGRATION_METHOD == NUMAV_ANALYTIC
            // calculate triangle area
            std::array<std::array<Float,DIM>,3UL> triangle_coords;
            for (size_t eni = 0UL; eni != 3UL; ++eni) {
                const size_t ni = _sei_to_ni[sei][eni];
                triangle_coords[eni] = std::array<Float,DIM>({
                    _ni_to_coords[ni][0UL],
                    _ni_to_coords[ni][1UL],
                    _ni_to_coords[ni][2UL]
                });
            }
            const Float triangle_area = get_triangle_area(triangle_coords);
        #endif

        // damping matrix
        #if NUMAV_DAMP_INTEGRATION_METHOD == NUMAV_GAUSS_QUADRATURE
            constexpr std::array<std::array<Float,2UL>,NGP_DAMP<O>>
                GAUSS_POINTS_DAMP = GAUSS_POINTS_SFC<NGP_DAMP<O>>;
            for (size_t gpi = 0UL; gpi != NGP_DAMP<O>; ++gpi)
            {
                const Eigen::Matrix<Float,ENI_COUNT_SFC<O>,2UL> nabla_n =
                    shape_func_sfc_gradient<O>(
                        GAUSS_POINTS_DAMP[gpi][0UL], GAUSS_POINTS_DAMP[gpi][1UL]
                    ); // todo: try putting constexpr here
                const Eigen::Matrix<Float,2UL,2UL> jac_matrix =
                    coords_matrix * nabla_n;
                const Float det_jac = std::abs(jac_matrix.determinant());
                
                const Eigen::Matrix<Float,ENI_COUNT_SFC<O>,1UL> n =
                shape_func_sfc<O>(
                    GAUSS_POINTS_DAMP[gpi][0UL], GAUSS_POINTS_DAMP[gpi][1UL]
                );
                
                const Eigen::Matrix<Float,ENI_COUNT_SFC<O>,ENI_COUNT_SFC<O>>
                    nnt = n * n.transpose();
                
                // todo: multiply detj and w without creating another matrix
                const Eigen::Matrix<Float,ENI_COUNT_SFC<O>,ENI_COUNT_SFC<O>>
                    nnt_detj_w =
                        nnt * det_jac * GAUSS_WEIGHTS_SFC<NGP_DAMP<O>>[gpi];
                
                for (size_t nci = 0UL; nci != ENI_PAIRS.size(); ++nci) {
                    _ispgi_to_damp_fi_part[ispgi][fipi_damp[nci]] += 
                        nnt_detj_w(ENI_PAIRS[nci][0UL], ENI_PAIRS[nci][1UL]);
                }
            }
        #elif NUMAV_DAMP_INTEGRATION_METHOD == NUMAV_ANALYTIC
            const Eigen::Matrix<Float,ENI_COUNT_SFC<O>,ENI_COUNT_SFC<O>>
                damp_matrix_const_part = DAMP_MATRIX_CONST_PART<O>;
            
            const Eigen::Matrix<Float,ENI_COUNT_SFC<O>,ENI_COUNT_SFC<O>>
                damp_fi_part = damp_matrix_const_part * triangle_area;
            
            for (size_t nci = 0UL; nci != ENI_PAIRS.size(); ++nci) {
                _ispgi_to_damp_fi_part[ispgi][fipi_damp[nci]] +=
                    damp_fi_part(ENI_PAIRS[nci][0UL], ENI_PAIRS[nci][1UL]);
            }
        #else
            static_assert(false, "Invalid NUMAV_DAMP_INTEGRATION_METHOD.");
        #endif
        }
    ispgi_to_map_to_fipi.free();
}

template<ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_assemble_fi_part_for_sfc_velocity()
{
    // count the fipi for each ispgv
    fz::SafePtr<std::unordered_map<size_t,size_t>> ispgv_to_map_to_fipi(
        _ispgv_count
    );
    for (size_t vsei = 0UL; vsei != _vsei_count; ++vsei) {
        const size_t ispgv = _vsei_to_ispgv[vsei];
        const size_t sei = _vsei_to_sei[vsei];
        for (size_t eni = 0UL; eni != ENI_COUNT_SFC<O>; ++eni) {
            const size_t ni = _sei_to_ni[sei][eni];
            if (!(ispgv_to_map_to_fipi[ispgv].count(ni) > 0)) {
                const size_t fipi = ispgv_to_map_to_fipi[ispgv].size();
                ispgv_to_map_to_fipi[ispgv].insert({ni, fipi});
            }
        }
    }

    // allocate memory in the safe pointers
    _ispgv_to_forc_fi_part = fz::SafePtr<fz::SafePtr<Float>>(_ispgv_count);
    _ispgv_to_ptr_in_b = fz::SafePtr<fz::SafePtr<Cmplx*>>(_ispgv_count);
    for (size_t ispgv = 0UL; ispgv != _ispgv_count; ++ispgv)
    {
        const size_t size = ispgv_to_map_to_fipi[ispgv].size();
        _ispgv_to_forc_fi_part[ispgv] = fz::SafePtr<Float>(size);
        _ispgv_to_forc_fi_part[ispgv].fill(0_F);
        _ispgv_to_ptr_in_b[ispgv] = fz::SafePtr<Cmplx*>(size);
    }

    // assemble elemental force vectors
    std::array<size_t, ENI_COUNT_SFC<O>> fipi_forc;
    for (size_t vsei = 0UL; vsei != _vsei_count; ++vsei)
    {
        const size_t ispgv = _vsei_to_ispgv[vsei];
        const size_t sei = _vsei_to_sei[vsei];

        // Create _ispgv_to_ptr_in_b
        for (size_t eni = 0UL; eni != ENI_COUNT_SFC<O>; ++eni)
        {
            const size_t ni = _sei_to_ni[sei][eni];
            fipi_forc[eni] = ispgv_to_map_to_fipi[ispgv].at(ni);
            
            const size_t* const ni_ptr = std::lower_bound(
                _b_row_idx.begin(), _b_row_idx.end(), ni
            );
            const size_t ptrdiff = ni_ptr - _b_row_idx.begin();
            _ispgv_to_ptr_in_b[ispgv][fipi_forc[eni]] =
                _b_vals.begin() + ptrdiff;
        }

        #if NUMAV_FORC_INTEGRATION_METHOD == NUMAV_GAUSS_QUADRATURE
            // coordinates matrix
            std::array<Eigen::Vector3d, ENI_COUNT_SFC<O>> triangle_3d;
            for (size_t eni = 0UL; eni != ENI_COUNT_SFC<O>; ++eni) {
                const size_t node_idx = _sei_to_ni[sei][eni];
                triangle_3d[eni] = Eigen::Vector3d(
                    _ni_to_coords[node_idx][0UL],
                    _ni_to_coords[node_idx][1UL],
                    _ni_to_coords[node_idx][2UL]
                );
            }
            std::array<Eigen::Vector2d, ENI_COUNT_SFC<O>> triangle_2d =
                project_triangle_to_2d<O>(triangle_3d);
            Eigen::Matrix<Float,2UL,ENI_COUNT_SFC<O>> coords_matrix;
            for (size_t eni = 0UL; eni != ENI_COUNT_SFC<O>; ++eni) {
                coords_matrix(0UL,eni) = triangle_2d[eni](0UL);
                coords_matrix(1UL,eni) = triangle_2d[eni](1UL);
            }
        #elif NUMAV_FORC_INTEGRATION_METHOD == NUMAV_ANALYTIC
            // calculate triangle area
            std::array<std::array<Float,DIM>,3UL> triangle_coords;
            for (size_t eni = 0UL; eni != 3UL; ++eni) {
                const size_t ni = _sei_to_ni[sei][eni];
                triangle_coords[eni] = std::array<Float,DIM>({
                    _ni_to_coords[ni][0UL],
                    _ni_to_coords[ni][1UL],
                    _ni_to_coords[ni][2UL]
                });
            }
            const Float triangle_area = get_triangle_area(triangle_coords);
        #endif

        // elemental force vector
        #if NUMAV_FORC_INTEGRATION_METHOD == NUMAV_GAUSS_QUADRATURE
            constexpr std::array<std::array<Float,2UL>,NGP_FORC<O>>
                GAUSS_POINTS_FORC = GAUSS_POINTS_SFC<NGP_FORC<O>>;
            for (size_t gpi = 0UL; gpi != NGP_FORC<O>; ++gpi)
            {
                const Eigen::Matrix<Float, ENI_COUNT_SFC<O>, 2UL> nabla_n =
                    shape_func_sfc_gradient<O>(
                        GAUSS_POINTS_FORC[gpi][0UL], GAUSS_POINTS_FORC[gpi][1UL]
                    ); // todo: try putting constexpr here
                const Eigen::Matrix<Float,2UL,2UL> jac_matrix =
                    coords_matrix * nabla_n;
                const Float det_jac = jac_matrix.determinant();

                const Eigen::Matrix<Float,ENI_COUNT_SFC<O>,1UL> n =
                    shape_func_sfc<O>(
                        GAUSS_POINTS_FORC[gpi][0UL], GAUSS_POINTS_FORC[gpi][1UL]
                    );

                // todo: multiply detj and w without creating another matrix
                const Eigen::Matrix<Float,ENI_COUNT_SFC<O>,1UL> n_detj_w =
                    n * det_jac * GAUSS_WEIGHTS_SFC<NGP_FORC<O>>[gpi];
                
                for (size_t eni = 0UL; eni != ENI_COUNT_SFC<O>; ++eni) {
                    _ispgv_to_forc_fi_part[ispgv][fipi_forc[eni]] +=
                        n_detj_w(eni);
                }
            }
        #elif NUMAV_FORC_INTEGRATION_METHOD == NUMAV_ANALYTIC
            const Eigen::Matrix<Float,ENI_COUNT_SFC<O>,1UL>
                forc_vector_const_part = FORC_VECTOR_CONST_PART<O>;
            
            const Eigen::Matrix<Float,ENI_COUNT_SFC<O>,1UL>
                forc_fi_part = forc_vector_const_part * triangle_area;
            
            for (size_t eni = 0UL; eni != ENI_COUNT_SFC<O>; ++eni) {
                _ispgv_to_forc_fi_part[ispgv][fipi_forc[eni]] +=
                    forc_fi_part(eni);
            }
        #else
            static_assert(false, "Invalid NUMAV_FORC_INTEGRATION_METHOD.");
        #endif
    }
    ispgv_to_map_to_fipi.free();
}

template<ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_assemble_fi_part_for_point_velocity()
{
    _vpi_to_ptr_in_b = fz::SafePtr<Cmplx*>(_vpi_count);
    for (size_t vpi = 0UL; vpi != _vpi_count; ++vpi)
    {
        const size_t ni = _vpi_to_ni[vpi];
        const size_t* const ni_ptr = std::lower_bound(
            _b_row_idx.begin(), _b_row_idx.end(), ni
        );
        const size_t ptrdiff = ni_ptr - _b_row_idx.begin();
        _vpi_to_ptr_in_b[vpi] = _b_vals.begin() + ptrdiff;
    }
}

template<typename T>
std::unordered_map<std::vector<size_t>, std::vector<T>> find_set_intersections(
    const fz::SafePtr<fz::SafePtr<T>>& sets
) {
    // map each element to the indices of sets that contain it
    std::unordered_map<T, std::vector<size_t>> element_to_sets;
    for (size_t set_index = 0UL; set_index != sets.size(); ++set_index) {
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
    fz::SafePtr<fz::SafePtr<size_t>> sets(_ppi_count + _ispgp_count);
    for (size_t ppi = 0UL; ppi != _ppi_count; ++ppi) {
        sets[ppi] = fz::SafePtr<size_t>(1UL); // TODO: review this
        sets[ppi][0UL] = _ppi_to_ni[ppi];
    }
    for (size_t ispgp = 0UL; ispgp != _ispgp_count; ++ispgp) {
        std::set<size_t> unique_nodes; // TODO: review the use of set
        for (size_t psei = 0UL; psei != _psei_count; ++psei) {
            if (_psei_to_ispgp[psei] == ispgp){
                const size_t sei = _psei_to_sei[psei];
                for (size_t eni = 0UL; eni != ENI_COUNT_SFC<O>; ++eni) {
                    const size_t ni = _sei_to_ni[sei][eni];
                    unique_nodes.insert(ni);
                }
            }
        }
        sets[_ppi_count+ispgp] = fz::SafePtr<size_t>(unique_nodes.size());
        size_t i = 0UL;
        for (const auto& ni : unique_nodes) {
            sets[_ppi_count+ispgp][i] = ni;
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
    _apvi_count = intersections.size();
    _apvi_to_pressure = fz::SafePtr<FuncFloatToCmplx>(_apvi_count);
    size_t apvi = 0UL;
    for (auto& [set_indexes, ni_vector] : intersections) {
        auto average = [this, set_indexes = set_indexes](const Float& freq) {
            Cmplx sum = Cmplx(0_F, 0_F);
            for (const auto& set_index : set_indexes) {
                if (set_index < _ppi_count) {
                    const size_t& ppi = set_index;
                    sum += (_ppi_to_pressure[ppi])(freq);
                }
                else {
                    const size_t ispgp = set_index - _ppi_count;
                    sum += (_ispgp_to_pressure[ispgp])(freq);
                }
            }
            return sum / static_cast<Float>(set_indexes.size());
        };
        _apvi_to_pressure[apvi] = average;
        ++apvi;
    }

    // define _apvi_to_ptr_in_a and _apvi_to_ptr_in_b
    _apvi_to_ptr_in_a = fz::SafePtr<fz::SafePtr<Cmplx*>>(_apvi_count);
    _apvi_to_ptr_in_b = fz::SafePtr<fz::SafePtr<Cmplx*>>(_apvi_count);
    apvi = 0UL;
    for (auto& [set_indexes, ni_vector] : intersections) {
        _apvi_to_ptr_in_a[apvi] = fz::SafePtr<Cmplx*>(ni_vector.size());
        _apvi_to_ptr_in_b[apvi] = fz::SafePtr<Cmplx*>(ni_vector.size());
        size_t fipi = 0UL;
        for (const auto& ni : ni_vector) {
            #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
                // _apvi_to_ptr_in_a
                _apvi_to_ptr_in_a[apvi][fipi] = _a_diag.data() + ni;
                const size_t* const ni_ptr = std::lower_bound(
                    _b_row_idx.begin(), _b_row_idx.end(), ni
                );
                // _apvi_to_ptr_in_b
                const size_t ptrdiff_b = ni_ptr - _b_row_idx.begin();
                _apvi_to_ptr_in_b[apvi][fipi] = _b_vals.begin() + ptrdiff_b;
                ++fipi;
                continue;
            #else
                // _apvi_to_ptr_in_a
                const std::pair<size_t,size_t> pair = std::make_pair(ni, ni);
                const std::pair<size_t,size_t>* const pair_ptr =
                    std::lower_bound(
                        _ni_connections.begin(),
                        _ni_connections.end(),
                        pair,
                        compare_pair<size_t>
                    );
                const size_t ptrdiff_a = 
                    pair_ptr - _ni_connections.begin();
                _apvi_to_ptr_in_a[apvi][fipi] = _a_vals.begin() + ptrdiff_a;
                
                // _apvi_to_ptr_in_b
                const size_t* const ni_ptr = std::lower_bound(
                    _b_row_idx.begin(), _b_row_idx.end(), ni
                );
                const size_t ptrdiff_b = ni_ptr - _b_row_idx.begin();
                _apvi_to_ptr_in_b[apvi][fipi] = _b_vals.begin() + ptrdiff_b;
                
                ++fipi;
            #endif
        }
        ++apvi;
    }
}

template<ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_assemble_freq_independent_parts()
{   
    _allocate_a();
    _allocate_b();
    _allocate_x();
    #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
        _define_sparsity_pattern_using_internal_solver();
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
        _define_sparsity_pattern_using_eigen_solver();
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
        _define_sparsity_pattern_using_onemkl_solver();
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_MUMPS
        _define_sparsity_pattern_using_mumps_solver();
    #else
        static_assert(false, "Invalid NUMAV_SYSTEM_SOLVER.");
    #endif
    _assemble_fi_part_for_vol_elements();
    _assemble_fi_part_for_sfc_impedance();
    _assemble_fi_part_for_point_velocity();
    _assemble_fi_part_for_sfc_velocity();
    _assemble_fi_part_for_pressure();
}

} // namespace numav

NUMAV_INSTANTIATE_SIM_AC_FEM_FREQ_D3
