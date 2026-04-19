// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "ldlt-solver.hpp"

#include <complex>
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream> // todo: remove this

template<typename T>
LdltSolver<T>::LdltSolver() = default;

template<typename T>
LdltSolver<T>::~LdltSolver() {
    delete[] _nzli_to_sum_count;
    delete[] _nzli_to_d;
    delete[] _nzli_to_a;
    delete[] _l_sum_factors;
    delete[] _l_row_idxs;
    delete[] _l_vals;
    delete[] _d;
    delete[] _nz_count_in_row;
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

    for (size_t ri=0; ri!=_row_count; ++ri) {
        if (a_row_idxs[ri] > a_row_idxs[ri+1]) {
            std::cout << "a_row_idxs is not in ascending order.\n";
            throw;
        }
    }

    const size_t nz_count = a_row_idxs[_row_count];
    for (size_t nzi=0; nzi!=nz_count; ++nzi) {
        if (a_col_idxs[nzi] > _row_count - 1) {
            std::cout << "a_col_idxs[" << nzi << "] is out of bounds.\n";
            throw;
        }
    }

    for (size_t ri=0; ri!=_row_count; ++ri) {
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

template<typename T>
void LdltSolver<T>::analyse_sparsity_pattern(
    const T* const a_diag,
    const size_t* const a_row_idxs,
    const size_t* const a_col_idxs,
    const T* const a_vals,
    const T* const b,
    const size_t& row_count
) {
    _a_diag = a_diag;
    _a_vals = a_vals;
    _b = b;
    _row_count = row_count;

    _check_if_input_is_valid(a_row_idxs, a_col_idxs);

    // todo: check if the system is diagonal

    // fist loop:
    //   Count sizes to allocate the containers that will be used soon.
    _l_row_idxs = new size_t[_row_count + 1];
    _l_row_idxs[0] = 0;
    size_t l_sum_factors_count = 0;
    for (size_t ri=0; ri!=_row_count; ++ri)
    {
        _l_row_idxs[ri+1] = _l_row_idxs[ri];
        for (size_t ci=0; ci!=ri; ++ci)
        {
            bool was_element_already_added = false;
            auto add_element_to_l_if_its_new = [&]() {
                if (!was_element_already_added) {
                    ++_l_row_idxs[ri+1];
                    _l_col_idxs.emplace_back(ci);
                    was_element_already_added = true;
                }
            };

            // check if a is nz for (ri,ci)
            const size_t a_row_idx = a_row_idxs[ri];
            const size_t a_row_idx_next = a_row_idxs[ri+1];
            const size_t* binary_search_begin = a_col_idxs + a_row_idx;
            const size_t* binary_search_end = a_col_idxs + a_row_idx_next;
            const bool is_a_nz = std::binary_search(
                binary_search_begin, binary_search_end, ci
            );
            if(is_a_nz) {
                add_element_to_l_if_its_new();
            }

            // check if the sum is nz for (ri,ci)
            if (_l_col_idxs.size() == 0) { continue; }
            size_t idx_it_1 = _l_row_idxs[ri];
            size_t idx_it_2 = _l_row_idxs[ci];
            size_t idx_it_1_end = _l_row_idxs[ri+1];
            size_t idx_it_2_end = _l_row_idxs[ci+1];
            while(
                _l_col_idxs[idx_it_1] < ci &&
                _l_col_idxs[idx_it_2] < ci &&
                idx_it_1 != idx_it_1_end &&
                idx_it_2 != idx_it_2_end
            ) {
                if (_l_col_idxs[idx_it_1] < _l_col_idxs[idx_it_2]) {
                    ++idx_it_1;
                }
                else if (_l_col_idxs[idx_it_1] > _l_col_idxs[idx_it_2]) {
                    ++idx_it_2;
                }
                else { // (*it_l_col_idxs_1 == *it_l_col_idxs_2)
                    add_element_to_l_if_its_new();
                    l_sum_factors_count += 3;
                    ++idx_it_1;
                    ++idx_it_2;
                }
            }
        }
    }
    _l_col_idxs.shrink_to_fit();
    const size_t nzli_count = _l_col_idxs.size();
    _l_vals = new T[nzli_count];
    _d = new T[_row_count];
    _nzli_to_d = new T*[nzli_count];
    _nzli_to_a = new const T*[nzli_count];
    _nzli_to_sum_count = new size_t[nzli_count];
    _l_sum_factors = new T* [l_sum_factors_count];
    
    // second loop:
    //   Get all the pointers organized in a way to speed up the solution as
    //   much as possible.
    T** it_nzli_to_d = _nzli_to_d;
    const T** it_nzli_to_a = _nzli_to_a;
    size_t* it_nzli_to_sum_count = _nzli_to_sum_count;
    T** it_l_sum_factors = _l_sum_factors;
    for (size_t ri=0; ri!=_row_count; ++ri)
    {
        for (size_t ci=0; ci!=ri; ++ci)
        {
            // check if a is nz for (ri,ci)
            const size_t a_row_idx = a_row_idxs[ri];
            const size_t a_row_idx_next = a_row_idxs[ri+1];
            const size_t* binary_search_begin = a_col_idxs + a_row_idx;
            const size_t* binary_search_end = a_col_idxs + a_row_idx_next;
            const size_t* ptr_to_nz_a = std::lower_bound(
                binary_search_begin, binary_search_end, ci
            );
            bool is_a_nz;
            if(ptr_to_nz_a!=binary_search_end && *ptr_to_nz_a==ci) {
                is_a_nz = true;
                size_t nzci = ptr_to_nz_a - binary_search_begin;
                *it_nzli_to_a = _a_vals + a_row_idxs[ri] + nzci;
                ++it_nzli_to_a;
            }
            else {
                is_a_nz = false;
            }

            // check if the sum is nz for (ri,ci)
            size_t sum_count = 0;
            const size_t* it_l_col_idxs_1 =
                _l_col_idxs.data() + _l_row_idxs[ri];
            const size_t* it_l_col_idxs_2 =
                _l_col_idxs.data() + _l_row_idxs[ci];
            const size_t* it_l_col_idxs_end_1 =
                _l_col_idxs.data() + _l_row_idxs[ri+1];
            const size_t* it_l_col_idxs_end_2 =
                _l_col_idxs.data() + _l_row_idxs[ci+1];
            while(
                *it_l_col_idxs_1 < ci &&
                *it_l_col_idxs_2 < ci &&
                it_l_col_idxs_1 != it_l_col_idxs_end_1 &&
                it_l_col_idxs_2 != it_l_col_idxs_end_2
            ) {
                if (*it_l_col_idxs_1 < *it_l_col_idxs_2) {
                    ++it_l_col_idxs_1;
                }
                else if (*it_l_col_idxs_1 > *it_l_col_idxs_2) {
                    ++it_l_col_idxs_2;
                }
                else { // (*it_l_col_idxs_1 == *it_l_col_idxs_2)
                    ++sum_count;

                    const size_t nzli_l1 = it_l_col_idxs_1 - _l_col_idxs.data();
                    T* const l1_ptr = _l_vals + nzli_l1;
                    const size_t nzli_l2 = it_l_col_idxs_2 - _l_col_idxs.data();
                    T* const l2_ptr = _l_vals+ nzli_l2;
                    const size_t ci_in_sum = *it_l_col_idxs_1;
                    T* const d_ptr = _d + ci_in_sum;

                    *it_l_sum_factors = l1_ptr;
                    ++it_l_sum_factors;
                    *it_l_sum_factors = l2_ptr;
                    ++it_l_sum_factors;
                    *it_l_sum_factors = d_ptr;
                    ++it_l_sum_factors;

                    ++it_l_col_idxs_1;
                    ++it_l_col_idxs_2;
                }
            }
            bool is_sum_nz = (sum_count != 0) ? true : false;

            if (!is_a_nz && is_sum_nz) {
                *it_nzli_to_a = nullptr;
                ++it_nzli_to_a;
            }
            if (is_a_nz || is_sum_nz) {
                *it_nzli_to_d = _d + ci;
                ++it_nzli_to_d;
                *it_nzli_to_sum_count = sum_count;
                ++it_nzli_to_sum_count;
            }
        }
    }

    // calculate the non zero count of elements in each row of L
    _nz_count_in_row = new size_t[_row_count];
    for (size_t ri=0; ri!=_row_count; ++ri)
    {
        const size_t row_ptr = _l_row_idxs[ri];
        const size_t next_row_ptr = _l_row_idxs[ri+1];
        _nz_count_in_row[ri] = next_row_ptr - row_ptr;
    }
}

template<typename T>
void LdltSolver<T>::solve(T* const x)
{
    T** it_nzli_to_d = _nzli_to_d;
    const T** it_nzli_to_a = _nzli_to_a;
    size_t* it_nzli_to_sum_count = _nzli_to_sum_count;
    T** it_l_sum_factors = _l_sum_factors;
    T* it_l_vals = _l_vals;
    for (size_t ri=0; ri!=_row_count; ++ri)
    {
        // lower matrix
        for (size_t nz_ci=0; nz_ci!=_nz_count_in_row[ri]; ++nz_ci)
        {
            // TODO: benchmark pre-computing the division
            T d_denom = **it_nzli_to_d;
            ++it_nzli_to_d;
            
            const size_t si_count = *it_nzli_to_sum_count;
            ++it_nzli_to_sum_count;
            if (si_count == 0) {
                *it_l_vals = **it_nzli_to_a / d_denom;
                ++it_l_vals;
                ++it_nzli_to_a;
                continue;
            }
            T sum = 0;
            for (size_t si=0; si!=si_count; ++si) {
                T product = **it_l_sum_factors;
                ++it_l_sum_factors;
                product *= **it_l_sum_factors;
                ++it_l_sum_factors;
                product *= **it_l_sum_factors;
                ++it_l_sum_factors;
                sum += product;
            }
            if (*it_nzli_to_a != nullptr) {
                *it_l_vals = (**it_nzli_to_a - sum) / d_denom;
                ++it_l_vals;
                ++it_nzli_to_a;
            }
            else {
                ++it_nzli_to_a;
                *it_l_vals = - sum / d_denom;
                ++it_l_vals;
            }
        }

        // diagonal
        const size_t row_ptr = _l_row_idxs[ri];
        T sum = 0;
        for (size_t nz_ci=0; nz_ci!=_nz_count_in_row[ri]; ++nz_ci)
        {
            const size_t l_ci_in_sum = _l_col_idxs[row_ptr + nz_ci];
            // TODO: benchmark other squaring methods
            sum += _d[l_ci_in_sum] * std::pow(_l_vals[row_ptr + nz_ci], 2);
        }
        _d[ri] = _a_diag[ri] - sum;
    }

    // compute y
    T* const& y = x;
    x[0] = _b[0];
    for(size_t ri=1; ri!=_row_count; ++ri) {
        T sum = 0;
        const size_t l_row_idx_begin = _l_row_idxs[ri];
        const size_t l_row_idx_end = _l_row_idxs[ri+1];
        for(size_t nzli=l_row_idx_begin; nzli!=l_row_idx_end; ++nzli) {
            const size_t ci = _l_col_idxs[nzli];
            sum += _l_vals[nzli] * y[ci];
        }
        y[ri] = _b[ri] - sum;
    }

    // compute z
    T* const& z = y;
    for(size_t ri=0; ri!=_row_count; ++ri) {
        z[ri] = y[ri] / _d[ri];
    }

    // compute x
    // x[_row_count-1] = z[_row_count-1]; (this line is not needed because x=z)
    for(size_t ri=_row_count-2; ri!=std::numeric_limits<size_t>::max(); --ri) {
        T sum = 0;
        for(size_t ci=ri+1; ci!=_row_count; ++ci) {
            const size_t* const l_row_idx_begin = &_l_col_idxs[_l_row_idxs[ci]];
            const size_t* const l_row_idx_end = &_l_col_idxs[_l_row_idxs[ci+1]];
            const size_t* const lb_ptr = std::lower_bound(
                l_row_idx_begin, l_row_idx_end, ri
            );
            if (lb_ptr!=l_row_idx_end && *lb_ptr==ri) {
                const size_t nzi = lb_ptr - _l_col_idxs.data();
                sum += _l_vals[nzi] * x[ci];
            }
        }
        x[ri] = z[ri] - sum;
    }
}

// Explicit template instantiations
template class LdltSolver<float>;
template class LdltSolver<double>;
template class LdltSolver<std::complex<float>>;
template class LdltSolver<std::complex<double>>;
