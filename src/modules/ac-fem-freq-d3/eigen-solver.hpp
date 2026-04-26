// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "SafePtr.hpp"
#include "numav/aliases.hpp"

namespace numav {

void solve_using_eigen(
    const fz::SafePtr<Cmplx>&,
    const fz::SafePtr<std::pair<size_t,size_t>>&,
    const fz::SafePtr<Cmplx>&,
    const fz::SafePtr<size_t>&,
    const size_t,
    Cmplx* const
);

} // namespace numav
