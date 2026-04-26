// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/eigen-solver.hpp"

#include "Eigen/Eigen"

namespace numav {

void solve_using_eigen(
    const fz::SafePtr<Cmplx>& a_vals,
    const fz::SafePtr<std::pair<size_t,size_t>>& nnz_rowcol_idx_pairs,
    const fz::SafePtr<Cmplx>& b_vals,
    const fz::SafePtr<size_t>& b_row_idx,
    const size_t node_count,
    Cmplx* const x_out
) {
    using Triplet = typename Eigen::Triplet<Cmplx>;

    // a matrix
    fz::SafePtr<Triplet> triplets_a(2UL*a_vals.size() - node_count);
    Triplet* it_triplets_a = triplets_a.begin();
    for (size_t j = 0UL; j != a_vals.size(); ++j) {
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
    Eigen::SparseMatrix<Cmplx> a(node_count, node_count);
    a.setFromTriplets(triplets_a.begin(), triplets_a.end());
    triplets_a.free();

    // b vector
    fz::SafePtr<Triplet> triplets_b(node_count);
    Triplet* it_triplets_b = triplets_b.begin();
    for (size_t j = 0UL; j != b_vals.size(); ++j) {
        *it_triplets_b = Triplet(b_row_idx[j], 0UL, b_vals[j]);
        ++it_triplets_b;
    }
    Eigen::SparseMatrix<Cmplx> b(node_count, 1UL);
    b.setFromTriplets(triplets_b.begin(), triplets_b.end());
    triplets_b.free();

    Eigen::SparseLU<Eigen::SparseMatrix<Cmplx>> solver;
    solver.analyzePattern(a);
    solver.factorize(a);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Factorization failed. Matrix may be singular.\n";
    }
    const Eigen::VectorXcd x = solver.solve(b);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solving failed.\n";
    }
    for (size_t j = 0UL; j != node_count; ++j) {
        x_out[j] = x(j);
    }
}

} // namespace numav
