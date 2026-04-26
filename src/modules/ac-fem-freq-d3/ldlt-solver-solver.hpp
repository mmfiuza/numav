// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/aliases.hpp"
#include "SafePtr.hpp"
#include "ldlt-solver.hpp"

#include <utility>

namespace numav {

void define_ldlt_solver_sparsity_pattern(
    LdltSolver<Cmplx>&,
    const fz::SafePtr<Cmplx>&,
    const fz::SafePtr<std::pair<size_t,size_t>>&,
    const fz::SafePtr<Cmplx>&,
    fz::SafePtr<Cmplx>&,
    fz::SafePtr<Cmplx>&,
    const size_t
);

void solve_using_ldlt_solver(
    LdltSolver<Cmplx>&,
    const fz::SafePtr<Cmplx>&,
    const fz::SafePtr<size_t>&,
    fz::SafePtr<Cmplx>&
);

} // namespace numav
