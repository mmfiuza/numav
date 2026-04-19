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

    void analyse_sparsity_pattern(
        const T* const,
        const size_t* const,
        const size_t* const,
        const T* const,
        const T* const,
        const size_t&
    );
    void solve(T* const);

private:
    void _check_if_input_is_valid(const size_t* const, const size_t* const);

    const T* _a_diag;
    const T* _a_vals;
    const T* _b;
    size_t _row_count;

    size_t* _nzli_to_sum_count;
    T** _nzli_to_d;
    const T** _nzli_to_a;
    T** _l_sum_factors;
    size_t* _l_row_idxs;
    std::vector<size_t> _l_col_idxs;
    T* _l_vals;
    T* _d;
    size_t* _nz_count_in_row;
};
