// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/exception.hpp"
#include "modules/ac-fem-freq-d3/impl.hpp"

#if NUMAV_SYSTEM_SOLVER == NUMAV_MUMPS

namespace numav {    

#define ONE_TO_ZERO_BASED_INDEX(I) (I-1)

constexpr MUMPS_INT JOB_INIT = -1;
constexpr MUMPS_INT JOB_END = -2;
constexpr MUMPS_INT JOB_ANALYSE = 1;
constexpr MUMPS_INT JOB_FACTORIZE = 2;
constexpr MUMPS_INT JOB_SOLVE = 3;

void print_mumps_error(const MUMPS_INT error_id) {
    error("MUMPS error code: {}", error_id);
}

template <ElementOrder O>
void
SimulationAcFemFreqD3Tet<O>::Impl::_define_sparsity_pattern_using_mumps_solver()
{
    // generate _a_row_idx and _a_col_idx
    const size_t nzi_count = _ni_connections.size();
    _a_row_idx = fz::SafePtr<MUMPS_INT>(nzi_count);
    _a_col_idx = fz::SafePtr<MUMPS_INT>(nzi_count);
    for (size_t nzi = 0UL; nzi != _ni_connections.size(); ++nzi) {
        _a_row_idx[nzi] = static_cast<MUMPS_INT>(
            _ni_connections[nzi].first + 1UL
        );
        _a_col_idx[nzi] = static_cast<MUMPS_INT>(
            _ni_connections[nzi].second + 1UL
        );
    }

    // allocate the dense b vector
    _b_dense = fz::SafePtr<Cmplx>(_ni_count);

    // initialize MUMPS
    _solver.job = JOB_INIT;
    _solver.par = 1; // host is working
    _solver.sym = 1;
    _solver.comm_fortran = 0;
    zmumps_c(&_solver);
    _solver.icntl[ONE_TO_ZERO_BASED_INDEX(1)] = 0; // no outputs
    _solver.icntl[ONE_TO_ZERO_BASED_INDEX(2)] = 0; // no outputs
    _solver.icntl[ONE_TO_ZERO_BASED_INDEX(3)] = 0; // no outputs
    _solver.icntl[ONE_TO_ZERO_BASED_INDEX(4)] = 0; // no outputs
    _solver.icntl[ONE_TO_ZERO_BASED_INDEX(5)] = 0; // assembled matrix
    _solver.icntl[ONE_TO_ZERO_BASED_INDEX(7)] = 3; // Scotch
    _solver.icntl[ONE_TO_ZERO_BASED_INDEX(20)] = 0; // dense rhs
    _solver.icntl[ONE_TO_ZERO_BASED_INDEX(28)] = 1; // sequential computation
    _solver.cntl[ONE_TO_ZERO_BASED_INDEX(1)] = 0_F; // no pivoting

    // pass input to solver
    _solver.n = static_cast<MUMPS_INT>(_ni_count);
    _solver.nnz = static_cast<MUMPS_INT>(nzi_count);
    _solver.irn = _a_row_idx.data();
    _solver.jcn = _a_col_idx.data();
    _solver.a = reinterpret_cast<mumps_double_complex*>(_a_vals.data());
    _solver.rhs = reinterpret_cast<mumps_double_complex*>(_b_dense.data());

    // analyze
    _solver.job = JOB_ANALYSE;
    zmumps_c(&_solver);
    if (_solver.infog[ONE_TO_ZERO_BASED_INDEX(1)] < 0) {
        print_mumps_error(_solver.infog[ONE_TO_ZERO_BASED_INDEX(2)]);
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3Tet<O>::Impl::_solve_using_mumps_solver()
{
    // factorize
    _solver.job = JOB_FACTORIZE;
    zmumps_c(&_solver);
    if (_solver.infog[ONE_TO_ZERO_BASED_INDEX(1)] < 0) {
        print_mumps_error(_solver.infog[ONE_TO_ZERO_BASED_INDEX(2)]);
    }

    // dense b vector
    _b_dense.fill(Cmplx(0_F, 0_F));
    for (size_t nzi = 0UL; nzi != _b_vals.size(); ++nzi) {
        _b_dense[_b_row_idx[nzi]] = _b_vals[nzi];
    }

    // solve
    _solver.job = JOB_SOLVE;
    zmumps_c(&_solver);
    if (_solver.infog[ONE_TO_ZERO_BASED_INDEX(1)] < 0) {
        print_mumps_error(_solver.infog[ONE_TO_ZERO_BASED_INDEX(2)]);
    }

    // copy solution
    for (size_t ni = 0UL; ni != _ni_count; ++ni) {
        _x[ni] = _b_dense[ni];
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3Tet<O>::Impl::_terminate_mumps_solver()
{
    _solver.job = JOB_END;
    zmumps_c(&_solver);
    _a_row_idx.free();
    _a_col_idx.free();
}

} // namespace numav

NUMAV_INSTANTIATE_SIM_AC_FEM_FREQ_D3

#endif
