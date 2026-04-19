// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/ldlt-solver-solver.hpp"

#include <limits>

namespace numav {

void define_ldlt_solver_sparsity_pattern(
    LdltSolver<_cmplx_t>& solver,
    const fz::SafePtr<_cmplx_t>& a_diag,
    const fz::SafePtr<std::pair<size_t,size_t>>& nz_rowcol_idx_pairs,
    const fz::SafePtr<_cmplx_t>& a_vals,
    fz::SafePtr<_cmplx_t>& b_dense,
    const size_t& ni_count
) {
    // problem dimensions
    const size_t nzi_count = nz_rowcol_idx_pairs.size();
    
    // create a_col_idx and a_row_idx
    fz::SafePtr<size_t> a_row_idx(ni_count + 1);
    fz::SafePtr<size_t> a_col_idx(nzi_count);
    size_t* it_a_row_idx = a_row_idx.begin();
    size_t* it_a_col_idx = a_col_idx.begin();
    const std::pair<size_t,size_t>* it_nz_rowcol_idx_pairs =
        nz_rowcol_idx_pairs.begin();
    size_t current_ri = std::numeric_limits<size_t>::max();
    for (size_t nzi=0; nzi!=nzi_count; ++nzi) {
        *it_a_col_idx = it_nz_rowcol_idx_pairs->second;
        ++it_a_col_idx;
        if (it_nz_rowcol_idx_pairs->first != current_ri) {
            const size_t previous_ri = current_ri;
            current_ri = it_nz_rowcol_idx_pairs->first;
            const size_t ri_increment = current_ri - previous_ri;
            for (size_t i=0; i!=ri_increment; ++i) {
                *it_a_row_idx = nzi;
                ++it_a_row_idx;
            }
        }
        ++it_nz_rowcol_idx_pairs;
    }
    *it_a_row_idx = nzi_count;

    // allocate the null dense b vector
    b_dense = fz::SafePtr<_cmplx_t>(ni_count);
    b_dense.fill(0.0);
    
    // define the non-zero structure of the matrix
    solver.analyse_sparsity_pattern(
        a_diag.data(),
        a_row_idx.data(),
        a_col_idx.data(),
        a_vals.data(),
        b_dense.data(),
        ni_count
    );

    a_row_idx.free();
    a_col_idx.free();
}

void solve_using_ldlt_solver(
    LdltSolver<_cmplx_t>& solver,
    const fz::SafePtr<_cmplx_t>& b_vals,
    const fz::SafePtr<size_t>& b_row_idx,
    fz::SafePtr<_cmplx_t>& b_dense,
    _cmplx_t* const x
) {
    // dense b vector
    for (size_t i=0; i!=b_vals.size(); ++i) {
        b_dense[b_row_idx[i]] = b_vals[i];
    }
    
    solver.solve(x);
}

} // namespace numav
