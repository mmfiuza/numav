// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/utils.hpp"
#include "modules/ac-fem-freq-d3/impl.hpp"
#include <limits>

#if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL

namespace numav {

template <ElementOrder O>
void
SimulationAcFemFreqD3Tet<O>::Impl::_define_sparsity_pattern_using_internal_solver()
{
    // problem dimensions
    const uint64_t nzi_count = _ni_connections.size();
    
    // create a_col_idx and a_row_idx
    fz::SafePtr<uint64_t> a_row_idx(_ni_count + 1UL);
    fz::SafePtr<uint64_t> a_col_idx(nzi_count);
    uint64_t* it_a_row_idx = a_row_idx.begin();
    uint64_t* it_a_col_idx = a_col_idx.begin();
    auto it_ni_connections = _ni_connections.begin();
    uint64_t current_ri = std::numeric_limits<uint64_t>::max();
    for (uint64_t nzi = 0UL; nzi != nzi_count; ++nzi) {
        *it_a_col_idx = it_ni_connections->second;
        ++it_a_col_idx;
        if (it_ni_connections->first != current_ri) {
            const uint64_t previous_ri = current_ri;
            current_ri = it_ni_connections->first;
            const uint64_t ri_increment = current_ri - previous_ri;
            for (uint64_t i = 0UL; i != ri_increment; ++i) {
                *it_a_row_idx = nzi;
                ++it_a_row_idx;
            }
        }
        ++it_ni_connections;
    }
    *it_a_row_idx = nzi_count;

    // allocate the null dense b vector
    _b_dense = fz::SafePtr<Cmplx>(_ni_count);
    _b_dense.fill(Cmplx(0_F, 0_F));
    
    // define the non-zero structure of the matrix
    _solver.define_sparsity_pattern(
        _a_diag.data(),
        a_row_idx.data(),
        a_col_idx.data(),
        _a_vals.data(),
        _x.data(),
        _b_dense.data(),
        _ni_count
    );

    a_row_idx.free();
    a_col_idx.free();
}
template <ElementOrder O>
void SimulationAcFemFreqD3Tet<O>::Impl::_solve_using_internal_solver()
{
    // dense b vector
    for (uint64_t i = 0UL; i != _b_vals.size(); ++i) {
        _b_dense[_b_row_idx[i]] = _b_vals[i];
    }
    
    _solver.solve();
}

} // namespace numav

NUMAV_INSTANTIATE_SIM_AC_FEM_FREQ_D3

#endif
