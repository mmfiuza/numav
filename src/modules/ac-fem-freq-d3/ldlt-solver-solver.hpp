// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "common/aliases.hpp"
#include "SafePtr.hpp"
#include "ldlt-solver.hpp"

#include <utility>

namespace numav {

void define_ldlt_solver_sparsity_pattern(
    LdltSolver<_cmplx_t>&,
    const fz::SafePtr<_cmplx_t>&,
    const fz::SafePtr<std::pair<size_t,size_t>>&,
    const fz::SafePtr<_cmplx_t>&,
    fz::SafePtr<_cmplx_t>&,
    fz::SafePtr<_cmplx_t>&,
    const size_t
);

void solve_using_ldlt_solver(
    LdltSolver<_cmplx_t>&,
    const fz::SafePtr<_cmplx_t>&,
    const fz::SafePtr<size_t>&,
    fz::SafePtr<_cmplx_t>&
);

} // namespace numav
