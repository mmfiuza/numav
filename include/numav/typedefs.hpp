// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#define SAFE_PTR_DEBUG
#undef NDEBUG

#include <cstddef>
#include <complex>
#include <functional>

// define the llinear system solver
#define NUMAV_EIGEN 01
#define NUMAV_ONEMKL 02
#define NUMAV_SYSTEM_SOLVER NUMAV_ONEMKL

// define the storage order for the global matrix
#define NUMAV_UPPER_COL_MAJOR 11
#define NUMAV_UPPER_ROW_MAJOR 12
#if NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
    #define NUMAV_GLOBAL_MATRIX_STORAGE_ORDER NUMAV_UPPER_COL_MAJOR
#elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
    #define NUMAV_GLOBAL_MATRIX_STORAGE_ORDER NUMAV_UPPER_ROW_MAJOR
#else
    static_assert(false, "Invalid NUMAV_SYSTEM_SOLVER.");
#endif

namespace numav {

    using _idx_t = size_t;
    using _epg_t = size_t;
    using _ipg_t = size_t;

    using _cmplx_t = typename std::complex<double>;

    using _FuncRealToCmplx = typename std::function<_cmplx_t(const double&)>;

} // namespace numav
