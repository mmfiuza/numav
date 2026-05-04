// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/exception.hpp"
#include "modules/ac-fem-freq-d3/impl.hpp"

namespace numav {

#if NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN

template <ElementOrder O>
void
SimulationAcFemFreqD3<O>::Impl::_define_sparsity_pattern_using_eigen_solver()
{
    // rewrite A matrix in CSC form
    const size_t nz_count = _ni_connections.size();
    _a_row_idx = fz::SafePtr<Eigen::Index>(nz_count);
    _a_col_idx = fz::SafePtr<Eigen::Index>(_ni_count + 1UL);
    Eigen::Index* it_a_row_idx = _a_row_idx.begin();
    Eigen::Index* it_a_col_idx = _a_col_idx.begin();
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

    // create the A matrix
    _a_eigen.emplace(
        _ni_count,
        _ni_count,
        nz_count,
        _a_col_idx.data(),
        _a_row_idx.data(),
        _a_vals.data()
    );

    // create b vector
    _b_col_idx_signed[0UL] = static_cast<Eigen::Index>(0);
    _b_col_idx_signed[1UL] = static_cast<Eigen::Index>(_b_vals.size());
    _b_row_idx_signed = fz::SafePtr<Eigen::Index>(_b_row_idx.size());
    for (size_t i = 0UL; i != _b_row_idx.size(); ++i) {
        _b_row_idx_signed[i] = static_cast<Eigen::Index>(_b_row_idx[i]);
    }
    _b_eigen.emplace(
        _ni_count,
        1UL,
        _b_vals.size(),
        _b_col_idx_signed.data(),
        _b_row_idx_signed.data(),
        _b_vals.data()
    );

    // create x vector
    _x_eigen.emplace(
        _x.data(), _ni_count
    );

    // analyze A sparsity pattern
    _solver->isSymmetric(true); // TODO: explore more what this does
    _solver->setPivotThreshold(Float(0.0));
    _solver->analyzePattern(*_a_eigen);
    if (_solver->info() != Eigen::Success) {
        error("Eigen::analyzePattern failed.");
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_solve_using_eigen_solver()
{
    _solver->factorize(*_a_eigen);
    if (_solver->info() != Eigen::Success) {
        error("Eigen::factorize failed.");
    }
    *_x_eigen = _solver->solve(*_b_eigen);
    if (_solver->info() != Eigen::Success) {
        error("Eigen::solve failed.");
    }
}

#endif

// explicit instantiation declarations
INSTANTIATE_SIMULATION_CLASS

} // namespace numav
