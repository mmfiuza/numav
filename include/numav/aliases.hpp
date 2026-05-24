// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include <functional>
#include <complex>
#include <array>

namespace numav {

    using Float = double;
    using Cmplx = typename std::complex<Float>;
    using FuncFloatToCmplx = typename std::function<Cmplx(const Float)>;
    using Coord = typename std::array<Float,3UL>;

} // namespace numav
