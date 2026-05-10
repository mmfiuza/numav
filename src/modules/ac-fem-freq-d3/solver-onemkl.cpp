// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/exception.hpp"
#include "common/utils.hpp"
#include "modules/ac-fem-freq-d3/impl.hpp"

#if NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL

namespace numav {


void print_dss_error(const MKL_INT error_id) {
    error("oneMLK error code: {}", error_id);
}

template <ElementOrder O>
void
SimulationAcFemFreqD3<O>::Impl::_define_sparsity_pattern_using_onemkl_solver()
{
    // problem dimensions
    const MKL_INT node_count = static_cast<MKL_INT>(_ni_count);
    const MKL_INT nnz_count = static_cast<MKL_INT>(_ni_connections.size());
    
    // create a_col_idx and a_row_ptr
    fz::SafePtr<MKL_INT> a_col_idx(nnz_count);
    fz::SafePtr<MKL_INT> a_row_ptr(node_count + 1UL);
    MKL_INT* it_a_col_idx = a_col_idx.begin();
    MKL_INT* it_a_row_ptr = a_row_ptr.begin();
    auto it_ni_connections = _ni_connections.begin();
    size_t current_row = std::numeric_limits<size_t>::max();
    for (MKL_INT i = static_cast<MKL_INT>(0); i != nnz_count; ++i) {
        *it_a_col_idx = it_ni_connections->second;
        ++it_a_col_idx;
        if (it_ni_connections->first != current_row) {
            current_row = it_ni_connections->first;
            *it_a_row_ptr = i; 
            ++it_a_row_ptr;
        }
        ++it_ni_connections;
    }
    *it_a_row_ptr = nnz_count;
    ++it_a_row_ptr;
    
    // error code
    MKL_INT error_id;
    
    // initialize the solver
    error_id = dss_create(_dss_handle, NUMAV_MKL_OPTIONS);
    if (error_id != MKL_DSS_SUCCESS) {
        print_dss_error(error_id);
    }
    
    // define the non-zero structure of the matrix
    const MKL_INT symmetry_type = MKL_DSS_SYMMETRIC_COMPLEX;
    error_id = dss_define_structure(
        _dss_handle,
        symmetry_type,
        a_row_ptr.data(),
        node_count,
        node_count,
        a_col_idx.data(),
        nnz_count
    );
    if (error_id != MKL_DSS_SUCCESS) {
        print_dss_error(error_id);
    }
    
    // reorder the matrix
    error_id = dss_reorder(_dss_handle, NUMAV_MKL_OPTIONS, 0UL);
    if (error_id != MKL_DSS_SUCCESS) {
        print_dss_error(error_id);
    }
    a_row_ptr.free();
    a_col_idx.free();

    // allocate the null dense b vector
    _b_dense = fz::SafePtr<Cmplx>(_ni_count);
    _b_dense.fill(Cmplx(0_F, 0_F));
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_solve_using_onemkl_solver()
{
    // error code
    MKL_INT error_id;
    
    // factor the matrix
    constexpr MKL_INT positive_definiteness = MKL_DSS_INDEFINITE;
    error_id = dss_factor_complex(
        _dss_handle,
        positive_definiteness,
        reinterpret_cast<const Float*>(_a_vals.data())
    );
    if (error_id != MKL_DSS_SUCCESS) {
        print_dss_error(error_id);
    }

    // dense b vector
    for (size_t i = 0UL; i != _b_vals.size(); ++i) {
        _b_dense[_b_row_idx[i]] = _b_vals[i];
    }
    
    // solve
    constexpr MKL_INT num_of_b = 1UL;
    error_id = dss_solve_complex(
        _dss_handle,
        NUMAV_MKL_OPTIONS,
        reinterpret_cast<const Float*>(_b_dense.data()),
        num_of_b,
        reinterpret_cast<Float*>(_x.data())
    );
    if (error_id != MKL_DSS_SUCCESS) {
        print_dss_error(error_id);
    }
}

// explicit instantiation declarations
INSTANTIATE_SIMULATION_CLASS

} // namespace numav

#endif
