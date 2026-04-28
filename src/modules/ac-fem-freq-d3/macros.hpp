// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/aliases.hpp"

#if NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
    #include "mkl_dss.h"
    #include "mkl_types.h"
#endif

// define the linear system solver
#define NUMAV_LDLT_SOLVER 11UL
#define NUMAV_EIGEN 12UL
#define NUMAV_ONEMKL 13UL
#define NUMAV_SYSTEM_SOLVER NUMAV_ONEMKL

// define the storage order for the global matrix
#define NUMAV_LOWER_ROW_MAJOR 21UL
#define NUMAV_LOWER_COL_MAJOR 22UL
#define NUMAV_UPPER_ROW_MAJOR 23UL
#define NUMAV_UPPER_COL_MAJOR 24UL
#if NUMAV_SYSTEM_SOLVER == NUMAV_LDLT_SOLVER
    #define NUMAV_GLOBAL_MATRIX_STORAGE_ORDER NUMAV_LOWER_ROW_MAJOR
#elif NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
    #define NUMAV_GLOBAL_MATRIX_STORAGE_ORDER NUMAV_UPPER_COL_MAJOR
#elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
    #define NUMAV_GLOBAL_MATRIX_STORAGE_ORDER NUMAV_UPPER_ROW_MAJOR
#else
    static_assert(false, "Invalid NUMAV_SYSTEM_SOLVER.");
#endif

// define the integration method for the matrices
#define NUMAV_ANALYTIC 31UL
#define NUMAV_GAUSS_QUADRATURE 32UL
#define NUMAV_STIF_INTEGRATION_METHOD NUMAV_GAUSS_QUADRATURE
#define NUMAV_MASS_INTEGRATION_METHOD NUMAV_ANALYTIC
#define NUMAV_DAMP_INTEGRATION_METHOD NUMAV_ANALYTIC
#define NUMAV_FORC_INTEGRATION_METHOD NUMAV_ANALYTIC

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

// explicit instantiation declarations
#define INSTANTIATE_SIMULATION_CLASS \
    template class Simulation< \
        Phenomenon::ACOUSTIC, \
        NumericalMethod::FEM, \
        Domain::FREQUENCY, \
        Dimension::D3, \
        ElementOrder::O1 \
    >; \
    template class Simulation< \
        Phenomenon::ACOUSTIC, \
        NumericalMethod::FEM, \
        Domain::FREQUENCY, \
        Dimension::D3, \
        ElementOrder::O2 \
    >;
