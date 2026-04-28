// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include <complex>
#include <functional>

namespace numav {

    using Float = float;
    using Cmplx = typename std::complex<Float>;
    using FuncFloatToCmplx = typename std::function<Cmplx(const double&)>;

} // namespace numav
