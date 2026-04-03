// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "common/aliases.hpp"
#include "common/debug-macros.hpp"

#if NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
    #include "mkl_dss.h"
    #include "mkl_types.h"
#endif

#include "SafePtr.hpp"

namespace numav {

#if NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
    void print_dss_error(const _INTEGER_t&);

    void define_onemkl_sparsity_pattern(
        _MKL_DSS_HANDLE_t&,
        const fz::SafePtr<std::pair<size_t,size_t>>&,
        const size_t,
        fz::SafePtr<_cmplx_t>&
    );

    void solve_using_onemkl(
        _MKL_DSS_HANDLE_t&,
        const fz::SafePtr<_cmplx_t>&,
        const fz::SafePtr<_cmplx_t>&,
        const fz::SafePtr<size_t>&,
        fz::SafePtr<_cmplx_t>&,
        _cmplx_t* const
    );

#endif

} // namespace numav
