// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/exception.hpp"
#include "common/utils.hpp"
#include "modules/ac-fem-freq-d3/onemkl-solver.hpp"
#include "modules/ac-fem-freq-d3/macros.hpp"

namespace numav {

#if NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL

    void print_dss_error(const _INTEGER_t& error_id) {
        error("oneMLK error code: {}", error_id);
    }

    void define_onemkl_sparsity_pattern(
        _MKL_DSS_HANDLE_t& dss_handle,
        const fz::SafePtr<std::pair<size_t,size_t>>& nnz_rowcol_idx_pairs,
        const size_t ni_count,
        fz::SafePtr<Cmplx>& b_dense
    ) {
        // problem dimensions
        const MKL_INT node_count = ni_count;
        const MKL_INT nnz_count = nnz_rowcol_idx_pairs.size();
        
        // create a_col_idx and a_row_ptr
        fz::SafePtr<MKL_INT> a_col_idx(nnz_count);
        fz::SafePtr<MKL_INT> a_row_ptr(node_count + 1UL);
        MKL_INT* it_a_col_idx = a_col_idx.begin();
        MKL_INT* it_a_row_ptr = a_row_ptr.begin();
        const std::pair<size_t,size_t>* it_nnz_rowcol_idx_pairs =
            nnz_rowcol_idx_pairs.begin();
        size_t current_row = std::numeric_limits<size_t>::max();
        for (MKL_INT i = 0; i != nnz_count; ++i) { // TODO: define type of 0
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
        
        // error code
        _INTEGER_t error_id;
        
        // initialize the solver
        error_id = dss_create(dss_handle, NUMAV_MKL_OPTIONS);
        if (error_id != MKL_DSS_SUCCESS) { print_dss_error(error_id); }
        
        // define the non-zero structure of the matrix
        const MKL_INT symmetry_type = MKL_DSS_SYMMETRIC_COMPLEX;
        error_id = dss_define_structure(
            dss_handle,
            symmetry_type,
            a_row_ptr.data(),
            node_count,
            node_count,
            a_col_idx.data(),
            nnz_count
        );
        if (error_id != MKL_DSS_SUCCESS) { print_dss_error(error_id); }
        
        // reorder the matrix
        error_id = dss_reorder(dss_handle, NUMAV_MKL_OPTIONS, 0UL);
        if (error_id != MKL_DSS_SUCCESS) { print_dss_error(error_id); }
        a_row_ptr.free();
        a_col_idx.free();

        // allocate the null dense b vector
        b_dense = fz::SafePtr<Cmplx>(ni_count);
        b_dense.fill(Cmplx(0_F, 0_F));
    }

    void solve_using_onemkl(
        _MKL_DSS_HANDLE_t& dss_handle,
        const fz::SafePtr<Cmplx>& a_vals,
        const fz::SafePtr<Cmplx>& b_vals,
        const fz::SafePtr<size_t>& b_row_idx,
        fz::SafePtr<Cmplx>& b_dense,
        Cmplx* const x_out
    ) {
        // error code
        _INTEGER_t error_id;
        
        // factor the matrix
        constexpr MKL_INT positive_definiteness = MKL_DSS_INDEFINITE;
        error_id = dss_factor_complex(
            dss_handle,
            positive_definiteness,
            reinterpret_cast<const Float*>(a_vals.data())
        );
        if (error_id != MKL_DSS_SUCCESS) { print_dss_error(error_id); }

        // dense b vector
        for (size_t i = 0UL; i != b_vals.size(); ++i) {
            b_dense[b_row_idx[i]] = b_vals[i];
        }
        
        // solve
        constexpr MKL_INT options = NUMAV_MKL_OPTIONS;
        constexpr MKL_INT num_of_b = 1UL;
        error_id = dss_solve_complex(
            dss_handle,
            options,
            reinterpret_cast<const Float*>(b_dense.data()),
            num_of_b,
            reinterpret_cast<Float*>(x_out)
        );
        if (error_id != MKL_DSS_SUCCESS) { print_dss_error(error_id); }
    }

#endif

} // namespace numav
