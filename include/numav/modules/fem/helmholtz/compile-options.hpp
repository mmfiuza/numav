// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

// define the linear system solver
#define NUMAV_INTERNAL 11
#define NUMAV_EIGEN 12
#define NUMAV_ONEMKL 13
#define NUMAV_MUMPS 14
#ifdef NUMAV_USE_INTERNAL_SOLVER
    #define NUMAV_SYSTEM_SOLVER NUMAV_INTERNAL
#elif defined(NUMAV_USE_EIGEN_SOLVER)
    #define NUMAV_SYSTEM_SOLVER NUMAV_EIGEN
#elif defined(NUMAV_USE_ONEMKL_SOLVER)
    #define NUMAV_SYSTEM_SOLVER NUMAV_ONEMKL
#elif defined(NUMAV_USE_MUMPS_SOLVER)
    #define NUMAV_SYSTEM_SOLVER NUMAV_MUMPS
#else
    static_assert(false, "Invalid NUMAV_SYSTEM_SOLVER.");
#endif

// include solvers if needed
#if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
    #include "ldlt-solver.hpp"
#elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
    #include "mkl_dss.h"
    #include "mkl_types.h"
#elif NUMAV_SYSTEM_SOLVER == NUMAV_MUMPS
    #include "zmumps_c.h"
#endif

// define the integration method for the matrices
#define NUMAV_ANALYTIC 21
#define NUMAV_GAUSS_QUADRATURE 22
#define NUMAV_STIF_INTEGRATION_METHOD NUMAV_ANALYTIC
#define NUMAV_MASS_INTEGRATION_METHOD NUMAV_ANALYTIC
#define NUMAV_DAMP_INTEGRATION_METHOD NUMAV_ANALYTIC
#define NUMAV_FORC_INTEGRATION_METHOD NUMAV_ANALYTIC

namespace numav {

// define the triangular type for the global matrix
enum class TriangularMatrixType : uint64_t {
    LOWER,
    UPPER
};
#if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
    constexpr TriangularMatrixType GLOBAL_MATRIX_TRIANGULAR_TYPE =
        TriangularMatrixType::LOWER;
#elif NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
    constexpr TriangularMatrixType GLOBAL_MATRIX_TRIANGULAR_TYPE =
        TriangularMatrixType::UPPER;
#elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
    constexpr TriangularMatrixType GLOBAL_MATRIX_TRIANGULAR_TYPE =
        TriangularMatrixType::UPPER;
#elif NUMAV_SYSTEM_SOLVER == NUMAV_MUMPS
    constexpr TriangularMatrixType GLOBAL_MATRIX_TRIANGULAR_TYPE =
        TriangularMatrixType::UPPER;
#else
    static_assert(false, "Invalid NUMAV_SYSTEM_SOLVER.");
#endif

// define the storage order for the global matrix
enum class MatrixStorageOrder : uint64_t {
    ROW_MAJOR,
    COL_MAJOR
};
#if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
    constexpr MatrixStorageOrder GLOBAL_MATRIX_STORAGE_ORDER =
        MatrixStorageOrder::ROW_MAJOR;
#elif NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
    constexpr MatrixStorageOrder GLOBAL_MATRIX_STORAGE_ORDER =
        MatrixStorageOrder::COL_MAJOR;
#elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
    constexpr MatrixStorageOrder GLOBAL_MATRIX_STORAGE_ORDER =
        MatrixStorageOrder::ROW_MAJOR;
#elif NUMAV_SYSTEM_SOLVER == NUMAV_MUMPS
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
