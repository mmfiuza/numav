// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/exception.hpp"
#include "modules/ac-fem-freq-d3/impl.hpp"

namespace numav {

#if NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN

template <ElementOrder O>
void
SimulationAcFemFreqD3<O>::Impl::_define_sparsity_pattern_using_eigen_solver()
{
    const size_t nz_count = _ni_connections.size();

    // rewrite matrix in CSC form
    _a_row_idx = fz::SafePtr<ptrdiff_t>(nz_count);
    _a_col_idx = fz::SafePtr<ptrdiff_t>(_ni_count + 1UL);
    ptrdiff_t* it_a_row_idx = _a_row_idx.begin();
    ptrdiff_t* it_a_col_idx = _a_col_idx.begin();
    auto it_ni_connections = _ni_connections.begin();
    size_t current_col = std::numeric_limits<size_t>::max();
    for (size_t nzi = 0UL; nzi != nz_count; ++nzi) {
        *it_a_row_idx = it_ni_connections->first;
        ++it_a_row_idx;
        if (it_ni_connections->second != current_col) {
            current_col = it_ni_connections->second;
            *it_a_col_idx = nzi; 
            ++it_a_col_idx;
        }
        ++it_ni_connections;
    }
    *it_a_col_idx = nz_count;

    _a.emplace(
        _ni_count,
        _ni_count,
        nz_count,
        _a_col_idx.data(),
        _a_row_idx.data(),
        _a_vals.data()
    );

    _solver->analyzePattern(*_a);
    if (_solver->info() != Eigen::Success) {
        error("Eigen::analyzePattern failed.");
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_solve_using_eigen_solver()
{
    using Triplet = typename Eigen::Triplet<Cmplx>;

    // b vector
    fz::SafePtr<Triplet> triplets_b(_ni_count);
    Triplet* it_triplets_b = triplets_b.begin();
    for (size_t j = 0UL; j != _b_vals.size(); ++j) {
        *it_triplets_b = Triplet(_b_row_idx[j], 0UL, _b_vals[j]);
        ++it_triplets_b;
    }
    Eigen::SparseMatrix<Cmplx> b(_ni_count, 1UL);
    b.setFromTriplets(triplets_b.begin(), triplets_b.end());
    triplets_b.free();

    _solver->factorize(*_a);
    if (_solver->info() != Eigen::Success) {
        error("Eigen::factorize failed.");
    }
    const Eigen::Matrix<Cmplx, Eigen::Dynamic, 1UL> x_temp = _solver->solve(b);
    if (_solver->info() != Eigen::Success) {
        error("Eigen::solve failed.");
    }
    for (size_t ni = 0UL; ni != _ni_count; ++ni) {
        _x[ni] = x_temp(ni);
    }
}

#endif

// explicit instantiation declarations
INSTANTIATE_SIMULATION_CLASS

} // namespace numav
