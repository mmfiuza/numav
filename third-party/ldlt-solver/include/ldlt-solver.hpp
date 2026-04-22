// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include <cstddef>
#include <memory>

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
    class Impl;
    std::unique_ptr<Impl> _pimpl;
};
