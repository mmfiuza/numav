// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include <cstddef>
#include <vector>

template<typename T>
class LdltSolver {
public:
    LdltSolver();
    ~LdltSolver();
    LdltSolver(const LdltSolver&) = delete;
    LdltSolver& operator=(const LdltSolver&) = delete;
    LdltSolver(LdltSolver&&) noexcept;
    LdltSolver& operator=(LdltSolver&&) noexcept;

    void define_sparsity_pattern(
        const T* const a_diag,
        const size_t* const a_row_idxs,
        const size_t* const a_col_idxs,
        const T* const a_vals,
        T* const x,
        const T* const b,
        const size_t& row_count
    );
    void solve();

private:
    void _check_if_input_is_valid(const size_t* const, const size_t* const);

    const T* __restrict__ _a_diag;
    const T* __restrict__ _a_vals;
    T* __restrict__ _x;
    const T* __restrict__ _b;
    size_t _unknown_count;
    size_t _li_count;
    size_t _nzi_count;

    T* __restrict__ _over_d;
    T* __restrict__ _l;
    T** __restrict__ _nzi_to_l_ptr;
};
