// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/numav.hpp"
#include "common/utils.hpp"
#include "modules/ac-fem-freq-d3/constants.hpp"
#include "Eigen/Eigen"

namespace numav {

template<ElementOrder O>
Eigen::Matrix<Float, ENIS_COUNT<O>, 1UL> shape_func_sfc(
    const Float xi0,
    const Float xi1
);

template<ElementOrder O>
Eigen::Matrix<Float, ENIS_COUNT<O>, 2UL> shape_func_sfc_gradient(
    const Float xi0,
    const Float xi1
);

template<ElementOrder O>
Eigen::Matrix<Float, ENIV_COUNT<O>, 1UL> shape_func_vol(
    const Float xi0,
    const Float xi1,
    const Float xi2
);

template<ElementOrder O>
Eigen::Matrix<Float, ENIV_COUNT<O>, DIM> shape_func_vol_gradient(
    const Float xi0,
    const Float xi1,
    const Float xi2
);

} // namespace numav
