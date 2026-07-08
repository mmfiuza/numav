// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "ldlt-solver.hpp"

#include <complex>
#include <cmath>
#include <algorithm>
#include <array>
#include <limits>
#include <iostream>

constexpr size_t TS = 32UL;

template<typename T>
LdltSolver<T>::LdltSolver() = default;

template<typename T>
LdltSolver<T>::~LdltSolver() {
    delete[] _over_d;
    delete[] _l;
    delete[] _nzi_to_l_ptr;
};

template<typename T>
LdltSolver<T>::LdltSolver(LdltSolver&& other) noexcept = default;

template<typename T>
LdltSolver<T>& LdltSolver<T>::operator=(LdltSolver&& other) noexcept = default;

template<typename T>
void LdltSolver<T>::_check_if_input_is_valid(
    const size_t* const a_row_idxs,
    const size_t* const a_col_idxs
) {
    if (a_row_idxs[0] != 0) {
        std::cout << "a_row_idxs first element is not zero.\n";
        throw;
    }

    for (size_t ri=0; ri!=_unknown_count; ++ri) {
        if (a_row_idxs[ri] > a_row_idxs[ri+1]) {
            std::cout << "a_row_idxs is not in ascending order.\n";
            throw;
        }
    }

    const size_t nz_count = a_row_idxs[_unknown_count];
    for (size_t nzi=0; nzi!=nz_count; ++nzi) {
        if (a_col_idxs[nzi] > _unknown_count - 1) {
            std::cout << "a_col_idxs[" << nzi << "] is out of bounds.\n";
            throw;
        }
    }

    for (size_t ri=0; ri!=_unknown_count; ++ri) {
        const size_t nzi_row_begin = a_row_idxs[ri];
        const size_t nzi_row_end = a_row_idxs[ri+1];
        if (nzi_row_begin == nzi_row_end) { continue; }
        for (size_t nzi=nzi_row_begin; nzi!=nzi_row_end-1; ++nzi) {
            if (a_col_idxs[nzi] >= a_col_idxs[nzi+1]) {
                std::cout << "a_col_idxs is not in ascending order.\n";
                throw;
            }
        }
    }
}

inline constexpr size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();

template<typename T>
void LdltSolver<T>::define_sparsity_pattern(
    const T* const a_diag,
    const size_t* const a_row_idx,
    const size_t* const a_col_idx,
    const T* const a_vals,
    T* const x,
    const T* const b,
    const size_t& unknown_count
) {
    _a_diag = a_diag;
    _a_vals = a_vals;
    _x = x;
    _b = b;
    _unknown_count = unknown_count;

    _check_if_input_is_valid(a_row_idx, a_col_idx);

    _li_count = _unknown_count * _unknown_count;
    _nzi_count = a_row_idx[_unknown_count];

    _l = new T[_li_count];
    _nzi_to_l_ptr = new T*[_nzi_count];
    
    const size_t* it_a_row_idx = a_row_idx;
    size_t nz_count_in_row = 0UL;
    size_t nzci = SIZE_T_MAX;
    size_t ri = SIZE_T_MAX;
    for (size_t nzi = 0UL; nzi != _nzi_count; ++nzi)
    {
        ++nzci;
        while (nzci == nz_count_in_row) {
            ++ri;
            nzci = 0UL;
            nz_count_in_row = *(it_a_row_idx + 1) - *it_a_row_idx;
            ++it_a_row_idx;
        }
        const size_t ci = a_col_idx[nzi];
        const size_t row_offset = ri;
        const size_t col_offset = ci*_unknown_count;
        const size_t offset = row_offset + col_offset;
        _nzi_to_l_ptr[nzi] = _l + offset;
    }

    _over_d = new T[_unknown_count];
}

template<typename T>
void LdltSolver<T>::solve()
{
    // initialize L as zero
    for (size_t li = 0UL; li != _li_count; ++li) {
        _l[li] = T(0);
    }
    
    // add A to L
    for (size_t nzi = 0UL; nzi != _nzi_count; ++nzi) {
        *_nzi_to_l_ptr[nzi] = _a_vals[nzi];
    }
    
    // add A to D
    for (size_t ri = 0UL; ri != _unknown_count; ++ri) {
        _over_d[ri] = _a_diag[ri];
    }
    
    const size_t& n = _unknown_count;

    // 1. Add to LD
    // 2. Calculate L
    // 3. Add to D
    // 4. Invert D
    for (size_t tci = 0UL; tci < n; tci += TS)
    {
        //! tri1 and tri2 in the first block
        for (size_t ci = tci; ci < std::min(tci+TS-1UL, n-1UL); ++ci)
        {
            const size_t ci_stride = ci*n;
            _over_d[ci] = T(1) / _over_d[ci]; //* invert D
            for (size_t ri1 = ci+1UL; ri1 < std::min(tci+TS-1UL, n-1UL); ++ri1)
            {
                const T l = _l[ri1 + ci_stride] * _over_d[ci]; //* calculate L
                _over_d[ri1] -= l * _l[ri1 + ci_stride]; //* add to L*LD to D
                _l[ri1 + ci_stride] = l; //* Update LD to L
                const size_t ri1_stride = ri1*n;
                for (size_t ri2 = ri1+1UL; ri2 < std::min(tci+TS, n); ++ri2)
                {
                    _l[ri2 + ri1_stride] -= _l[ri1 + ci_stride] * _l[ri2 + ci_stride]; //* add to L*LD to LD
                }
            }
            //* remaining ri1
            const size_t rem_ri1 = std::min(tci+TS-1UL, n-1UL);
            const T l = _l[rem_ri1 + ci_stride] * _over_d[ci]; //* calculate L
            _over_d[rem_ri1] -= l * _l[rem_ri1 + ci_stride]; //* add to L*LD to D
            _l[rem_ri1 + ci_stride] = l; //* Update LD to L
        }
        //* remaining ci
        const size_t rem_ci = std::min(tci+TS-1UL, n-1UL);
        _over_d[rem_ci] = T(1) / _over_d[rem_ci]; //* invert last D of block
        //! tri1 in the first block and tri2 in all others
        for (size_t tri2 = tci+TS; tri2 < n; tri2 += TS)
        {
            for (size_t ci = tci; ci < std::min(tci+TS-1UL, n-1UL); ++ci)
            {
                const size_t ci_stride = ci*n;
                for (size_t ri1 = ci+1UL; ri1 < std::min(tci+TS, n); ++ri1)
                {
                    const size_t ri1_stride = ri1*n;
                    for (size_t ri2 = tri2; ri2 < std::min(tri2+TS, n); ++ri2)
                    {
                        _l[ri2 + ri1_stride] -= _l[ri1 + ci_stride] * _l[ri2 + ci_stride]; //* add to L*LD to LD
                    }
                }
            }
        }
        for (size_t tri1 = tci+TS; tri1 < n; tri1 += TS)
        {
            //! every first iteration of tri1 in a given tci
            for (size_t ci = tci; ci < std::min(tci+TS, n-1UL); ++ci)
            {
                const size_t ci_stride = ci*n;
                for (size_t ri1 = tri1; ri1 < std::min(tri1+TS-1UL, n-1UL); ++ri1)
                {
                    const T l = _l[ri1 + ci_stride] * _over_d[ci]; //* calculate L
                    _over_d[ri1] -= l * _l[ri1 + ci_stride]; //* add to L*LD to D
                    _l[ri1 + ci_stride] = l; //* Update LD to L
                    const size_t ri1_stride = ri1*n;
                    for (size_t ri2 = ri1+1UL; ri2 < std::min(tri1+TS, n); ++ri2)
                    {
                        _l[ri2 + ri1_stride] -= _l[ri1 + ci_stride] * _l[ri2 + ci_stride]; //* add to L*LD to LD
                    }
                }
                //* remaining ri1
                const size_t rem_ri1 = std::min(tri1+TS-1UL, n-1UL);
                const T l = _l[rem_ri1 + ci_stride] * _over_d[ci]; //* calculate L
                _over_d[rem_ri1] -= l * _l[rem_ri1 + ci_stride]; //* add to L*LD to D
                _l[rem_ri1 + ci_stride] = l; //* Update LD to L
            }

            //! square block of tri1 and tri2
            for (size_t tri2 = tri1+TS; tri2 < n; tri2 += TS)
            {
                for (size_t ci = tci; ci < std::min(tci+TS, n-1UL); ++ci)
                {
                    const size_t ci_stride = ci*n;
                    for (size_t ri1 = tri1; ri1 < std::min(tri1+TS, n-1UL); ++ri1)
                    {
                        const size_t ri1_stride = ri1*n;
                        for (size_t ri2 = tri2; ri2 < std::min(tri2+TS, n); ++ri2)
                        {
                            _l[ri2 + ri1_stride] -= _l[ri1 + ci_stride] * _l[ri2 + ci_stride]; //* add to L*LD to LD
                        }
                    }
                }
            }
        }
    }

    // compute y
    for(size_t ui = 0UL; ui != n; ++ui) {
        _x[ui] = _b[ui];
    }
    for(size_t ci = 0UL; ci != n-1UL; ++ci)
    {
        const size_t ci_stride = ci*n;
        for(size_t ri = ci+1UL; ri != n; ++ri)
        {
            _x[ri] -= _l[ri + ci_stride] * _x[ci];
        }
    }

    // compute z
    for(size_t ui = 0; ui != n; ++ui) {
        _x[ui] *= _over_d[ui];
    }

    // compute x
    // x[_row_count-1] = z[_row_count-1]; (this line is not needed because x=z)
    for(size_t ui = n-2UL; ui != SIZE_T_MAX; --ui) {
        const size_t ui_stride = ui*n;
        for(size_t si = ui+1UL; si != n; ++si) {
            _x[ui] -= _l[si + ui_stride] * _x[si];
        }
    }
}

// Explicit template instantiations
template class LdltSolver<float>;
template class LdltSolver<double>;
template class LdltSolver<std::complex<float>>;
template class LdltSolver<std::complex<double>>;
