// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

// define the linear system solver
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

// define the integration method for the matrices
#define NUMAV_ANALYTIC 31
#define NUMAV_GAUSS_QUADRATURE 32
#define NUMAV_STIF_INTEGRATION_METHOD NUMAV_GAUSS_QUADRATURE
#define NUMAV_MASS_INTEGRATION_METHOD NUMAV_ANALYTIC
#define NUMAV_DAMP_INTEGRATION_METHOD NUMAV_ANALYTIC

#define NUMAV_MKL_OPTIONS \
    MKL_DSS_MSG_LVL_WARNING + \
    MKL_DSS_TERM_LVL_ERROR + \
    MKL_DSS_ZERO_BASED_INDEXING
