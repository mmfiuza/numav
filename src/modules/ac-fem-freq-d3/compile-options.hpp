// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/aliases.hpp"

// define the linear system solver
#define NUMAV_LDLT_SOLVER 11
#define NUMAV_EIGEN 12
#define NUMAV_ONEMKL 13
#define NUMAV_SYSTEM_SOLVER NUMAV_ONEMKL

// include solvers if needed
#if NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
    #include "mkl_dss.h"
    #include "mkl_types.h"
#endif

// define the integration method for the matrices
#define NUMAV_ANALYTIC 21
#define NUMAV_GAUSS_QUADRATURE 22
#define NUMAV_STIF_INTEGRATION_METHOD NUMAV_GAUSS_QUADRATURE
#define NUMAV_MASS_INTEGRATION_METHOD NUMAV_ANALYTIC
#define NUMAV_DAMP_INTEGRATION_METHOD NUMAV_ANALYTIC
#define NUMAV_FORC_INTEGRATION_METHOD NUMAV_ANALYTIC

namespace numav {

// define the triangular type for the global matrix
enum class TriangularMatrixType {
    LOWER,
    UPPER
};
#if NUMAV_SYSTEM_SOLVER == NUMAV_LDLT_SOLVER
    constexpr TriangularMatrixType GLOBAL_MATRIX_TRIANGULAR_TYPE =
        TriangularMatrixType::LOWER;
#elif NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
    constexpr TriangularMatrixType GLOBAL_MATRIX_TRIANGULAR_TYPE =
        TriangularMatrixType::UPPER;
#elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
    constexpr TriangularMatrixType GLOBAL_MATRIX_TRIANGULAR_TYPE =
        TriangularMatrixType::UPPER;
#else
    static_assert(false, "Invalid NUMAV_SYSTEM_SOLVER.");
#endif

// define the storage order for the global matrix
enum class MatrixStorageOrder {
    ROW_MAJOR,
    COL_MAJOR
};
#if NUMAV_SYSTEM_SOLVER == NUMAV_LDLT_SOLVER
    constexpr MatrixStorageOrder GLOBAL_MATRIX_STORAGE_ORDER =
        MatrixStorageOrder::ROW_MAJOR;
#elif NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
    constexpr MatrixStorageOrder GLOBAL_MATRIX_STORAGE_ORDER =
        MatrixStorageOrder::COL_MAJOR;
#elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
    constexpr MatrixStorageOrder GLOBAL_MATRIX_STORAGE_ORDER =
        MatrixStorageOrder::ROW_MAJOR;
#else
    static_assert(false, "Invalid NUMAV_SYSTEM_SOLVER.");
#endif

// Intel oneMKL solver options
#if NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
    constexpr MKL_INT NUMAV_MKL_OPTIONS =
        []<typename F = numav::Float>() constexpr {
            MKL_INT options =
                MKL_DSS_MSG_LVL_WARNING +
                MKL_DSS_TERM_LVL_ERROR +
                MKL_DSS_ZERO_BASED_INDEXING;
            
            if constexpr (std::is_same_v<F, float>) {
                options += MKL_DSS_SINGLE_PRECISION;
            }
            return options;
        }();
#endif

} // namespace numav
