// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#define SAFE_PTR_DEBUG
#undef NDEBUG

#include <cstddef>
#include <complex>
#include <functional>

namespace numav {

    using _idx_t = size_t;
    using _epg_t = size_t;
    using _ipg_t = size_t;

    using _cmplx_t = typename std::complex<double>;

    using _FuncRealToCmplx = typename std::function<_cmplx_t(const double&)>;

} // namespace numav
