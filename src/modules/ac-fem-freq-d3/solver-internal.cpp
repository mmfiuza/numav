// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/utils.hpp"

#include "modules/ac-fem-freq-d3/impl.hpp"

#include <limits>

namespace numav {

#if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL

template <ElementOrder O>
void
SimulationAcFemFreqD3<O>::Impl::define_sparsity_pattern_using_internal_solver()
{
    // problem dimensions
    const size_t nzi_count = _ni_connections.size();
    
    // create a_col_idx and a_row_idx
    fz::SafePtr<size_t> a_row_idx(_ni_count + 1UL);
    fz::SafePtr<size_t> a_col_idx(nzi_count);
    size_t* it_a_row_idx = a_row_idx.begin();
    size_t* it_a_col_idx = a_col_idx.begin();
    auto it_ni_connections = _ni_connections.begin();
    size_t current_ri = std::numeric_limits<size_t>::max();
    for (size_t nzi = 0UL; nzi != nzi_count; ++nzi) {
        *it_a_col_idx = it_ni_connections->second;
        ++it_a_col_idx;
        if (it_ni_connections->first != current_ri) {
            const size_t previous_ri = current_ri;
            current_ri = it_ni_connections->first;
            const size_t ri_increment = current_ri - previous_ri;
            for (size_t i = 0UL; i != ri_increment; ++i) {
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
void SimulationAcFemFreqD3<O>::Impl::solve_using_internal_solver()
{
    // dense b vector
    for (size_t i = 0UL; i != _b_vals.size(); ++i) {
        _b_dense[_b_row_idx[i]] = _b_vals[i];
    }
    
    _solver.solve();
}

#endif

// explicit instantiation declarations
INSTANTIATE_SIMULATION_CLASS

} // namespace numav
