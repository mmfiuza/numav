// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include <complex>
#include <functional>

namespace numav {

    using _cmplx_t = typename std::complex<double>;
    using _FuncRealToCmplx = typename std::function<_cmplx_t(const double&)>;

} // namespace numav
