// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "ldlt-solver.hpp"

#include <complex>
#include <algorithm>
#include <iostream>
#include <vector>
#include <limits>

template<typename T>
class LdltSolver<T>::Impl
{
public:
    Impl() = default;
    ~Impl();

    void define_sparsity_pattern(
        const T* const a_diag,
        const size_t* const a_row_idxs,
        const size_t* const a_col_idxs,
        const T* const a_vals,
        T* const x,
        const T* const b,
        const size_t& row_count
    );
    void solve();

private:
    void _error(
        const char* const message
    );
    void _check_if_input_is_valid(
        const size_t* const a_row_idxs,
        const size_t* const a_col_idxs
    );
    void _check_if_a_is_diagonal(
        const size_t* const a_row_idxs
    );

    struct _NzilData {
        size_t sum_count;
        T* over_d;
        const T* a;
    };

    const T* _a_diag;
    const T* _a_vals;
    T* _x;
    const T* _b;
    size_t _row_count;
    const T** _l_sum_factors;
    size_t* _l_row_idx;
    std::vector<size_t> _l_col_idx;
    T* _l_vals;
    T* _ld_vals;
    T* _over_d;
    size_t* _nzil_count_in_row;
    size_t* _nzil_count_in_col;
    _NzilData* _nzil_data;
    std::vector<const T*> _y_sum_factors;
    std::vector<const T*> _x_sum_factors;
    bool _is_sparsity_pattern_defined = false;
    bool _is_a_diagonal;
};

template<typename T>
LdltSolver<T>::Impl::~Impl() {
    if (_is_sparsity_pattern_defined && !_is_a_diagonal) {
        delete[] _l_sum_factors;
        delete[] _l_row_idx;
        delete[] _l_vals;
        delete[] _ld_vals;
        delete[] _over_d;
        delete[] _nzil_count_in_row;
        delete[] _nzil_count_in_col;
        delete[] _nzil_data;
    }
}

template<typename T>
void LdltSolver<T>::Impl::_error(
    const char* const message
) {
    std::cerr <<
    "[" << "\033[31m" << "LDLT Solver error" << "\033[0m" << "]: " <<
    message << "\n";
    throw std::runtime_error("LDLT Solver error");
}

template<typename T>
void LdltSolver<T>::Impl::_check_if_input_is_valid(
    const size_t* const a_row_idxs,
    const size_t* const a_col_idxs
) {
    if (a_row_idxs[0] != 0) {
        _error("a_row_idx's first element is not zero.");
    }

    for (size_t ri = 0; ri != _row_count; ++ri) {
        if (a_row_idxs[ri] > a_row_idxs[ri + 1]) {
            _error("a_row_idx is not in ascending order.");
        }
    }

    const size_t nz_count = a_row_idxs[_row_count];
    for (size_t nzi = 0; nzi != nz_count; ++nzi) {
        if (a_col_idxs[nzi] > _row_count - 1) {
            _error("Some element in a_col_idx is out of bounds.");
        }
    }

    for (size_t ri = 0; ri != _row_count; ++ri) {
        const size_t nzi_row_begin = a_row_idxs[ri];
        const size_t nzi_row_end = a_row_idxs[ri + 1];
        if (nzi_row_begin == nzi_row_end) { continue; }
        for (size_t nzi = nzi_row_begin; nzi != nzi_row_end - 1; ++nzi) {
            if (a_col_idxs[nzi] >= a_col_idxs[nzi + 1]) {
                _error("a_col_idx is not in ascending order for every row.");
            }
        }
    }
}

template<typename T>
void LdltSolver<T>::Impl::_check_if_a_is_diagonal(
    const size_t* const a_row_idxs
) {
    _is_a_diagonal = true;
    for (size_t ri = 0; ri != _row_count; ++ri) {
        if (a_row_idxs[ri] != 0) {
            _is_a_diagonal = false;
            return;
        }
    }
}

inline constexpr size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();

template<typename T>
void LdltSolver<T>::Impl::define_sparsity_pattern(
    const T* const a_diag,
    const size_t* const a_row_idxs,
    const size_t* const a_col_idxs,
    const T* const a_vals,
    T* const x,
    const T* const b,
    const size_t& row_count
) {
    _a_diag = a_diag;
    _a_vals = a_vals;
    _x = x;
    _b = b;
    _row_count = row_count;

    if (_is_sparsity_pattern_defined) {
        _error("Sparsity pattern is already defined.");
    }

    _check_if_input_is_valid(a_row_idxs, a_col_idxs);

    _check_if_a_is_diagonal(a_row_idxs);
    if (_is_a_diagonal) {
        _is_sparsity_pattern_defined = true;
        return;
    }

    // fist loop:
    //   Count sizes to allocate the containers that will be used soon.
    _l_row_idx = new size_t[_row_count + 1];
    _l_row_idx[0] = 0;
    size_t l_sum_factors_count = 0;
    for (size_t ri = 0; ri != _row_count; ++ri) {
        _l_row_idx[ri + 1] = _l_row_idx[ri];
        for (size_t ci = 0; ci != ri; ++ci) {
            bool was_element_already_added = false;
            auto add_element_to_l_if_its_new = [&]() {
                if (!was_element_already_added) {
                    ++_l_row_idx[ri + 1];
                    _l_col_idx.emplace_back(ci);
                    was_element_already_added = true;
                }
            };

            // check if a is nz for (ri,ci)
            const size_t a_row_idx = a_row_idxs[ri];
            const size_t a_row_idx_next = a_row_idxs[ri + 1];
            const size_t* binary_search_begin = a_col_idxs + a_row_idx;
            const size_t* binary_search_end = a_col_idxs + a_row_idx_next;
            const bool is_a_nz = std::binary_search(
                binary_search_begin, binary_search_end, ci
            );
            if (is_a_nz) {
                add_element_to_l_if_its_new();
            }

            // check if the sum is nz for (ri,ci)
            if (_l_col_idx.size() == 0) { continue; }
            size_t idx_it_1 = _l_row_idx[ri];
            size_t idx_it_2 = _l_row_idx[ci];
            size_t idx_it_1_end = _l_row_idx[ri + 1];
            size_t idx_it_2_end = _l_row_idx[ci + 1];
            while (
                _l_col_idx[idx_it_1] < ci &&
                _l_col_idx[idx_it_2] < ci &&
                idx_it_1 != idx_it_1_end &&
                idx_it_2 != idx_it_2_end
            ) {
                if (_l_col_idx[idx_it_1] < _l_col_idx[idx_it_2]) {
                    ++idx_it_1;
                } else if (_l_col_idx[idx_it_1] > _l_col_idx[idx_it_2]) {
                    ++idx_it_2;
                } else { // (*it_l_col_idxs_1 == *it_l_col_idxs_2)
                    add_element_to_l_if_its_new();
                    l_sum_factors_count += 2;
                    ++idx_it_1;
                    ++idx_it_2;
                }
            }
        }
    }
    _l_col_idx.shrink_to_fit();
    const size_t nzil_count = _l_col_idx.size();
    _nzil_data = new _NzilData[nzil_count];
    _l_vals = new T[nzil_count];
    _ld_vals = new T[nzil_count];
    _over_d = new T[_row_count];
    _l_sum_factors = new const T*[l_sum_factors_count];

    // second loop:
    //   Get all the pointers organized in a way to speed up the solution as
    //   much as possible.
    _NzilData* it_nzil_data = _nzil_data;
    const T** it_l_sum_factors = _l_sum_factors;
    for (size_t ri = 0; ri != _row_count; ++ri) {
        for (size_t ci = 0; ci != ri; ++ci) {
            // check if a is nz for (ri,ci)
            const size_t a_row_idx = a_row_idxs[ri];
            const size_t a_row_idx_next = a_row_idxs[ri + 1];
            const size_t* binary_search_begin = a_col_idxs + a_row_idx;
            const size_t* binary_search_end = a_col_idxs + a_row_idx_next;
            const size_t* ptr_to_nz_a = std::lower_bound(
                binary_search_begin, binary_search_end, ci
            );
            bool is_a_nz;
            if (ptr_to_nz_a != binary_search_end && *ptr_to_nz_a == ci) {
                is_a_nz = true;
                size_t nzci = ptr_to_nz_a - binary_search_begin;
                it_nzil_data->a = _a_vals + a_row_idxs[ri] + nzci;
            } else {
                is_a_nz = false;
            }

            // check if the sum is nz for (ri,ci)
            size_t sum_count = 0;
            const size_t* it_l_col_idxs_1 =
                _l_col_idx.data() + _l_row_idx[ri];
            const size_t* it_l_col_idxs_2 =
                _l_col_idx.data() + _l_row_idx[ci];
            const size_t* it_l_col_idxs_end_1 =
                _l_col_idx.data() + _l_row_idx[ri + 1];
            const size_t* it_l_col_idxs_end_2 =
                _l_col_idx.data() + _l_row_idx[ci + 1];
            while (
                *it_l_col_idxs_1 < ci &&
                *it_l_col_idxs_2 < ci &&
                it_l_col_idxs_1 != it_l_col_idxs_end_1 &&
                it_l_col_idxs_2 != it_l_col_idxs_end_2
            ) {
                if (*it_l_col_idxs_1 < *it_l_col_idxs_2) {
                    ++it_l_col_idxs_1;
                } else if (*it_l_col_idxs_1 > *it_l_col_idxs_2) {
                    ++it_l_col_idxs_2;
                } else { // (*it_l_col_idxs_1 == *it_l_col_idxs_2)
                    ++sum_count;

                    const size_t nzil_l = it_l_col_idxs_1 - _l_col_idx.data();
                    T* const l_ptr = _l_vals + nzil_l;
                    const size_t nzil_ld = it_l_col_idxs_2 - _l_col_idx.data();
                    T* const ld_ptr = _ld_vals + nzil_ld;

                    *it_l_sum_factors = l_ptr;
                    ++it_l_sum_factors;
                    *it_l_sum_factors = ld_ptr;
                    ++it_l_sum_factors;

                    ++it_l_col_idxs_1;
                    ++it_l_col_idxs_2;
                }
            }
            bool is_sum_nz = (sum_count != 0) ? true : false;

            if (!is_a_nz && is_sum_nz) {
                it_nzil_data->a = nullptr;
            }
            if (is_a_nz || is_sum_nz) {
                it_nzil_data->over_d = _over_d + ci;
                it_nzil_data->sum_count = sum_count;
                ++it_nzil_data;
            }
        }
    }

    // calculate the non zero count of elements in each row of L
    _nzil_count_in_row = new size_t[_row_count];
    for (size_t ri = 0; ri != _row_count; ++ri) {
        const size_t row_idx = _l_row_idx[ri];
        const size_t next_row_idx = _l_row_idx[ri + 1];
        _nzil_count_in_row[ri] = next_row_idx - row_idx;
    }

    // pre-compute pointers to sum factors of y
    T* const& y = _x;
    for (size_t ri = 1; ri != _row_count; ++ri) {
        const size_t l_row_idx_begin = _l_row_idx[ri];
        const size_t l_row_idx_end = _l_row_idx[ri + 1];
        for (size_t nzi = l_row_idx_begin; nzi != l_row_idx_end; ++nzi) {
            _y_sum_factors.emplace_back(_l_vals + nzi);
            _y_sum_factors.emplace_back(y + _l_col_idx[nzi]);
        }
    }
    _y_sum_factors.shrink_to_fit();

    // pre-compute pointers to sum factors of x
    _nzil_count_in_col = new size_t[_row_count];
    for (size_t ci = 0; ci != _row_count; ++ci) {
        _nzil_count_in_col[ci] = 0;
    }
    for (size_t ci = _row_count - 2; ci != SIZE_T_MAX; --ci) {
        for (size_t ri = ci + 1; ri != _row_count; ++ri) {
            const size_t* const l_row_idx_begin = &_l_col_idx[_l_row_idx[ri]];
            const size_t* const l_row_idx_end = &_l_col_idx[_l_row_idx[ri + 1]];
            const size_t* const lb_ptr = std::lower_bound(
                l_row_idx_begin, l_row_idx_end, ci
            );
            if (lb_ptr != l_row_idx_end && *lb_ptr == ci) {
                ++_nzil_count_in_col[ci];
                const size_t nzi = lb_ptr - _l_col_idx.data();
                _x_sum_factors.emplace_back(_l_vals + nzi);
                _x_sum_factors.emplace_back(_x + ri);
            }
        }
    }
    _x_sum_factors.shrink_to_fit();

    _is_sparsity_pattern_defined = true;
}

template<typename T>
void LdltSolver<T>::Impl::solve()
{
    if (!_is_sparsity_pattern_defined) {
        _error("Sparsity pattern not defined.");
    }

    if (_is_a_diagonal) {
        for (size_t ri = 0; ri != _row_count; ++ri) {
            _x[ri] = _b[ri] / _a_diag[ri];
        }
        return;
    }

    _NzilData* it_nzil_data = _nzil_data;
    T* it_l_vals = _l_vals;
    T* it_ld_vals = _ld_vals;
    const T** it_l_sum_factors = _l_sum_factors;
    const T* it_d_sum_factors_l = _l_vals;
    const T* it_d_sum_factors_ld = _ld_vals;
    for (size_t ri = 0; ri != _row_count; ++ri) {
        // lower matrix
        for (size_t nzci = 0; nzci != _nzil_count_in_row[ri]; ++nzci) {
            const size_t si_count = it_nzil_data->sum_count;
            T sum = 0;
            for (size_t si = 0; si != si_count; ++si) {
                T product = **it_l_sum_factors;
                ++it_l_sum_factors;
                product *= **it_l_sum_factors;
                ++it_l_sum_factors;
                sum += product;
            }
            if (it_nzil_data->a != nullptr) {
                *it_ld_vals = *it_nzil_data->a - sum;
            } else {
                *it_ld_vals = -sum;
            }
            const T over_d = *it_nzil_data->over_d;
            *it_l_vals = *it_ld_vals * over_d;
            ++it_l_vals;
            ++it_ld_vals;
            ++it_nzil_data;
        }

        // diagonal matrix
        T sum = 0;
        for (size_t nzci = 0; nzci != _nzil_count_in_row[ri]; ++nzci) {
            T product = *it_d_sum_factors_l;
            ++it_d_sum_factors_l;
            product *= *it_d_sum_factors_ld;
            ++it_d_sum_factors_ld;
            sum += product;
        }
        _over_d[ri] = T(1.0) / (_a_diag[ri] - sum);
    }

    // compute y
    T* const& y = _x;
    y[0] = _b[0];
    auto it_y_sum_factors = _y_sum_factors.begin();
    for (size_t ri = 1; ri != _row_count; ++ri) {
        T sum = 0;
        for (size_t nzci = 0; nzci != _nzil_count_in_row[ri]; ++nzci) {
            T product = **it_y_sum_factors;
            ++it_y_sum_factors;
            product *= **it_y_sum_factors;
            ++it_y_sum_factors;
            sum += product;
        }
        y[ri] = _b[ri] - sum;
    }

    // compute z
    T* const& z = y;
    for (size_t ri = 0; ri != _row_count; ++ri) {
        z[ri] = y[ri] * _over_d[ri];
    }

    // compute x
    // x[_row_count-1] = z[_row_count-1]; (this line is not needed because x=z)
    auto it_x_sum_factors = _x_sum_factors.begin();
    for (size_t ri = _row_count - 2; ri != SIZE_T_MAX; --ri) {
        const size_t& ci_l = ri;
        T sum = 0;
        for (size_t nzci = 0; nzci != _nzil_count_in_col[ci_l]; ++nzci) {
            T product = **it_x_sum_factors;
            ++it_x_sum_factors;
            product *= **it_x_sum_factors;
            ++it_x_sum_factors;
            sum += product;
        }
        _x[ri] = z[ri] - sum;
    }
}

// public interface – forwards to _pimpl
template<typename T>
LdltSolver<T>::LdltSolver() {
    _pimpl = std::make_unique<Impl>();
}
template<typename T>
LdltSolver<T>::~LdltSolver() = default;
template<typename T>
LdltSolver<T>::LdltSolver(LdltSolver&&) noexcept = default;
template<typename T>
LdltSolver<T>& LdltSolver<T>::operator=(LdltSolver&&) noexcept = default;
template<typename T>
void LdltSolver<T>::define_sparsity_pattern(
    const T* const a_diag,
    const size_t* const a_row_idxs,
    const size_t* const a_col_idxs,
    const T* const a_vals,
    T* const x,
    const T* const b,
    const size_t& row_count
) {
    _pimpl->define_sparsity_pattern(
        a_diag,
        a_row_idxs,
        a_col_idxs,
        a_vals,
        x,
        b,
        row_count
    );
}
template<typename T>
void LdltSolver<T>::solve() {
    _pimpl->solve();
}

// explicit template instantiations
template class LdltSolver<float>;
template class LdltSolver<double>;
template class LdltSolver<std::complex<float>>;
template class LdltSolver<std::complex<double>>;
