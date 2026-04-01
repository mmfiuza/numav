// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/impl.hpp"

#include "common/exception.hpp"
#include "common/maths.hpp"
#include "common/hash-functions.hpp"
#include "common/utils.hpp"
#include "modules/ac-fem-freq-d3/shape-functions.hpp"
#include "modules/ac-fem-freq-d3/analytic-integration.hpp"
#include "modules/ac-fem-freq-d3/gauss-quadrature.hpp"
#include "modules/ac-fem-freq-d3/onemkl-solver.hpp"

#include <set>

namespace numav {

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
void SimulationAcFemFreqD3<O>::Impl::_allocate_a()
{
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
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_allocate_b()
{
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
        #if \
        NUMAV_STIF_INTEGRATION_METHOD == NUMAV_GAUSS_QUADRATURE || \
        NUMAV_MASS_INTEGRATION_METHOD == NUMAV_GAUSS_QUADRATURE
            Eigen::Matrix<double,NODES_IN_VOL_ELEM<O>,DIM> coords_matrix;
            for (size_t eni=0; eni!=NODES_IN_VOL_ELEM<O>; ++eni) {
                const size_t ni = _vei_to_ni[vei][eni];
                for (size_t di=0; di!=DIM; ++di) {
                    coords_matrix(eni,di) = _ni_to_coords[ni][di];
                }
            }
        #endif
        #if \
        NUMAV_STIF_INTEGRATION_METHOD == NUMAV_ANALYTIC || \
        NUMAV_MASS_INTEGRATION_METHOD == NUMAV_ANALYTIC
            std::array<std::array<double,3>,4> coords;
            for (size_t eni=0; eni!=4; ++eni) {
                const size_t ni = _vei_to_ni[vei][eni];
                for (size_t di=0; di!=DIM; ++di) {
                    coords[eni][di] = _ni_to_coords[ni][di];
                }
            }
            const double tet_volume = get_tetrahedron_volume(coords);
        #endif

        // stiffness matrix
        #if NUMAV_STIF_INTEGRATION_METHOD == NUMAV_GAUSS_QUADRATURE
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

                const Eigen::Matrix<double,DIM,DIM> inv_jac =
                    jac_matrix.inverse();
                
                const Eigen::Matrix<double,DIM,NODES_IN_VOL_ELEM<O>> b_matrix =
                    inv_jac * nabla_n;
                
                const
                Eigen::Matrix<double,NODES_IN_VOL_ELEM<O>,NODES_IN_VOL_ELEM<O>>
                    btb = b_matrix.transpose() * b_matrix;
                
                const double det_jac = std::abs(jac_matrix.determinant());
                
                // todo: multiply detj and w without creating another matrix
                const
                Eigen::Matrix<double,NODES_IN_VOL_ELEM<O>,NODES_IN_VOL_ELEM<O>>
                    btb_detj_w = 
                        btb * det_jac * GAUSS_WEIGHTS_VOL<NGP_STIF<O>>[gpi];

                for (size_t nci=0; nci!=COMBS_VOL.size(); ++nci) {
                    _ivpg_to_stif_fi_part[ivpg][fipi_vol[nci]] +=
                        btb_detj_w(COMBS_VOL[nci][0], COMBS_VOL[nci][1]);
                }
            }
        #elif NUMAV_STIF_INTEGRATION_METHOD == NUMAV_ANALYTIC
            const
            Eigen::Matrix<double, NODES_IN_VOL_ELEM<O>, NODES_IN_VOL_ELEM<O>>
                stif_matrix_const_part = get_stif_matrix_const_part<O>(coords);

            const
            Eigen::Matrix<double, NODES_IN_VOL_ELEM<O>, NODES_IN_VOL_ELEM<O>>
                stif_fi_part = stif_matrix_const_part / (36*tet_volume);

            for (size_t nci=0; nci!=COMBS_VOL.size(); ++nci) {
                _ivpg_to_stif_fi_part[ivpg][fipi_vol[nci]] +=
                    stif_fi_part(COMBS_VOL[nci][0], COMBS_VOL[nci][1]);
            }
        #endif

        // mass matrix
        #if NUMAV_MASS_INTEGRATION_METHOD == NUMAV_GAUSS_QUADRATURE
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
                
                const double det_jac = std::abs(jac_matrix.determinant());

                const Eigen::Matrix<double,NODES_IN_VOL_ELEM<O>,1> n =
                    shape_func_vol<O>(
                        GAUSS_POINTS_MASS[gpi][0],
                        GAUSS_POINTS_MASS[gpi][1],
                        GAUSS_POINTS_MASS[gpi][2]
                    );

                const
                Eigen::Matrix<double,NODES_IN_VOL_ELEM<O>,NODES_IN_VOL_ELEM<O>>
                    nnt = n * n.transpose();
                
                // todo: multiply detj and w without creating another matrix
                const
                Eigen::Matrix<double,NODES_IN_VOL_ELEM<O>,NODES_IN_VOL_ELEM<O>>
                    nnt_detj_w =
                        nnt * det_jac * GAUSS_WEIGHTS_VOL<NGP_MASS<O>>[gpi];
                
                for (size_t nci=0; nci!=COMBS_VOL.size(); ++nci) {
                    _ivpg_to_mass_fi_part[ivpg][fipi_vol[nci]] +=
                        nnt_detj_w(COMBS_VOL[nci][0], COMBS_VOL[nci][1]);
                }
            }
        #elif NUMAV_MASS_INTEGRATION_METHOD == NUMAV_ANALYTIC
            const 
            Eigen::Matrix<double,NODES_IN_VOL_ELEM<O>,NODES_IN_VOL_ELEM<O>>
                mass_matrix_const_part = MASS_MATRIX_CONST_PART<O>;
            
            const
            Eigen::Matrix<double,NODES_IN_VOL_ELEM<O>,NODES_IN_VOL_ELEM<O>>
                mass_fi_part = mass_matrix_const_part * tet_volume;
            
            for (size_t nci=0; nci!=COMBS_VOL.size(); ++nci) {
                _ivpg_to_mass_fi_part[ivpg][fipi_vol[nci]] +=
                    mass_fi_part(COMBS_VOL[nci][0], COMBS_VOL[nci][1]);
            }
        #endif
    }
    ivpg_to_map_to_fipi.free();
}

template<ElementOrder O>
std::array<Eigen::Vector2d, NODES_IN_SFC_ELEM<O>> project_triangle_to_2d(
    const std::array<Eigen::Vector3d, NODES_IN_SFC_ELEM<O>>& vertices_3d
) {
    std::array<Eigen::Vector3d, NODES_IN_SFC_ELEM<O>> point_minus_origin;
    for (size_t eni = 0; eni!=NODES_IN_SFC_ELEM<O>; ++eni) {
        point_minus_origin[eni] = vertices_3d[eni] - vertices_3d[0];
    }
    const Eigen::Vector3d& u = point_minus_origin[1];
    const Eigen::Vector3d& v = point_minus_origin[2];
    const Eigen::Vector3d n = u.cross(v);
    const Eigen::Vector3d x = u / u.norm();
    Eigen::Vector3d y = n.cross(u);
    y /= y.norm();

    std::array<Eigen::Vector2d, NODES_IN_SFC_ELEM<O>> vertices_2d;
    for (size_t eni = 0; eni!=NODES_IN_SFC_ELEM<O>; ++eni) {
        vertices_2d[eni] = Eigen::Vector2d(
            point_minus_origin[eni].dot(x), point_minus_origin[eni].dot(y)
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

        #if NUMAV_DAMP_INTEGRATION_METHOD == NUMAV_GAUSS_QUADRATURE
            // coordinates matrix
            std::array<Eigen::Vector3d, NODES_IN_SFC_ELEM<O>> triangle_3d;
            for (size_t eni=0; eni!=NODES_IN_SFC_ELEM<O>; ++eni) {
                const size_t ni = _sei_to_ni[sei][eni];
                triangle_3d[eni] = Eigen::Vector3d(
                    _ni_to_coords[ni][0],
                    _ni_to_coords[ni][1],
                    _ni_to_coords[ni][2]
                );
            }
            std::array<Eigen::Vector2d, NODES_IN_SFC_ELEM<O>> triangle_2d =
                project_triangle_to_2d<O>(triangle_3d);
            Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,2> coords_matrix;
            for (size_t eni=0; eni!=NODES_IN_SFC_ELEM<O>; ++eni) {
                coords_matrix(eni,0) = triangle_2d[eni](0);
                coords_matrix(eni,1) = triangle_2d[eni](1);
            }
        #elif NUMAV_DAMP_INTEGRATION_METHOD == NUMAV_ANALYTIC
            // calculate triangle area
            std::array<std::array<double,DIM>,3> triangle_coords;
            for (size_t eni=0; eni!=3; ++eni) {
                const size_t ni = _sei_to_ni[sei][eni];
                triangle_coords[eni] = std::array<double,DIM>({
                    _ni_to_coords[ni][0],
                    _ni_to_coords[ni][1],
                    _ni_to_coords[ni][2]
                });
            }
            const double triangle_area = get_triangle_area(triangle_coords);
        #else
            static_assert(false, "Invalid NUMAV_DAMP_INTEGRATION_METHOD.");
        #endif

        // damping matrix
        #if NUMAV_DAMP_INTEGRATION_METHOD == NUMAV_GAUSS_QUADRATURE
            constexpr std::array<std::array<double,2>,NGP_DAMP<O>>
                GAUSS_POINTS_DAMP = GAUSS_POINTS_SFC<NGP_DAMP<O>>;
            for (size_t gpi=0; gpi!=NGP_DAMP<O>; ++gpi)
            {
                const Eigen::Matrix<double,2,NODES_IN_SFC_ELEM<O>> nabla_n =
                    shape_func_sfc_gradient<O>(
                        GAUSS_POINTS_DAMP[gpi][0], GAUSS_POINTS_DAMP[gpi][1]
                    ); // todo: try putting constexpr here
                const Eigen::Matrix<double,2,2> jac_matrix =
                    nabla_n * coords_matrix;
                const double det_jac = std::abs(jac_matrix.determinant());
                
                const Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,1> n =
                shape_func_sfc<O>(
                    GAUSS_POINTS_DAMP[gpi][0], GAUSS_POINTS_DAMP[gpi][1]
                );
                
                const
                Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,NODES_IN_SFC_ELEM<O>>
                    nnt = n * n.transpose();
                
                // todo: multiply detj and w without creating another matrix
                const
                Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,NODES_IN_SFC_ELEM<O>>
                    nnt_detj_w =
                        nnt * det_jac * GAUSS_WEIGHTS_SFC<NGP_DAMP<O>>[gpi];
                
                for (size_t nci=0; nci!=COMBS_SFC.size(); ++nci) {
                    _ispgi_to_damp_fi_part[ispgi][fipi_damp[nci]] += 
                        nnt_detj_w(COMBS_SFC[nci][0], COMBS_SFC[nci][1]);
                }
            }
        #elif NUMAV_DAMP_INTEGRATION_METHOD == NUMAV_ANALYTIC
            const 
            Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,NODES_IN_SFC_ELEM<O>>
                damp_matrix_const_part = DAMP_MATRIX_CONST_PART<O>;
            
            const
            Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,NODES_IN_SFC_ELEM<O>>
                damp_fi_part = damp_matrix_const_part * triangle_area;
            
            for (size_t nci=0; nci!=COMBS_SFC.size(); ++nci) {
                _ispgi_to_damp_fi_part[ispgi][fipi_damp[nci]] +=
                    damp_fi_part(COMBS_SFC[nci][0], COMBS_SFC[nci][1]);
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

        #if NUMAV_FORC_INTEGRATION_METHOD == NUMAV_GAUSS_QUADRATURE
            // coordinates matrix
            std::array<Eigen::Vector3d, NODES_IN_SFC_ELEM<O>> triangle_3d;
            for (size_t eni=0; eni!=NODES_IN_SFC_ELEM<O>; ++eni) {
                const size_t node_idx = _sei_to_ni[sei][eni];
                triangle_3d[eni] = Eigen::Vector3d(
                    _ni_to_coords[node_idx][0],
                    _ni_to_coords[node_idx][1],
                    _ni_to_coords[node_idx][2]
                );
            }
            std::array<Eigen::Vector2d, NODES_IN_SFC_ELEM<O>> triangle_2d =
                project_triangle_to_2d<O>(triangle_3d);
            Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,2> coords_matrix;
            for (size_t eni=0; eni!=NODES_IN_SFC_ELEM<O>; ++eni) {
                coords_matrix(eni,0) = triangle_2d[eni](0);
                coords_matrix(eni,1) = triangle_2d[eni](1);
            }
        #elif NUMAV_FORC_INTEGRATION_METHOD == NUMAV_ANALYTIC
            // calculate triangle area
            std::array<std::array<double,DIM>,3> triangle_coords;
            for (size_t eni=0; eni!=3; ++eni) {
                const size_t ni = _sei_to_ni[sei][eni];
                triangle_coords[eni] = std::array<double,DIM>({
                    _ni_to_coords[ni][0],
                    _ni_to_coords[ni][1],
                    _ni_to_coords[ni][2]
                });
            }
            const double triangle_area = get_triangle_area(triangle_coords);
        #else
            static_assert(false, "Invalid NUMAV_FORC_INTEGRATION_METHOD.");
        #endif

        // elementary force vector
        #if NUMAV_FORC_INTEGRATION_METHOD == NUMAV_GAUSS_QUADRATURE
            constexpr std::array<std::array<double,2>,NGP_FORC<O>>
                GAUSS_POINTS_FORC = GAUSS_POINTS_SFC<NGP_FORC<O>>;
            for (size_t gpi=0; gpi!=NGP_FORC<O>; ++gpi)
            {
                const Eigen::Matrix<double,2,NODES_IN_SFC_ELEM<O>> nabla_n =
                    shape_func_sfc_gradient<O>(
                        GAUSS_POINTS_FORC[gpi][0], GAUSS_POINTS_FORC[gpi][1]
                    ); // todo: try putting constexpr here
                const Eigen::Matrix<double,2,2> jac_matrix =
                    nabla_n * coords_matrix;
                const double det_jac = jac_matrix.determinant();

                const Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,1> n =
                    shape_func_sfc<O>(
                        GAUSS_POINTS_FORC[gpi][0], GAUSS_POINTS_FORC[gpi][1]
                    );

                // todo: multiply detj and w without creating another matrix
                const Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,1> n_detj_w =
                    n * det_jac * GAUSS_WEIGHTS_SFC<NGP_FORC<O>>[gpi];
                
                for (size_t eni=0; eni!=NODES_IN_SFC_ELEM<O>; ++eni) {
                    _ispgv_to_forc_fi_part[ispgv][fipi_forc[eni]] +=
                        n_detj_w(eni);
                }
            }
        #elif NUMAV_FORC_INTEGRATION_METHOD == NUMAV_ANALYTIC
            const 
            Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,1>
                forc_vector_const_part = FORC_VECTOR_CONST_PART<O>;
            
            const
            Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,1>
                forc_fi_part = forc_vector_const_part * triangle_area;
            
            for (size_t eni=0; eni!=NODES_IN_SFC_ELEM<O>; ++eni) {
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
        std::set<size_t> unique_nodes; // TODO: review the use of set
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
    _allocate_a();
    _allocate_b();
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

// explicit instantiation declarations
INSTANTIATE_SIMULATION_CLASS

} // namespace numav
