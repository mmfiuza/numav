// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/numav.hpp"
#include "common/utils.hpp"
#include "modules/ac-fem-freq-d3/constants.hpp"
#include "Eigen/Eigen"

namespace numav {

template<ElementOrder O>
Eigen::Matrix<Float, NODES_IN_SFC_ELEM<O>, 1UL> shape_func_sfc(
    const Float,
    const Float
);

template<ElementOrder O>
Eigen::Matrix<Float, 2UL, NODES_IN_SFC_ELEM<O>> shape_func_sfc_gradient(
    const Float,
    const Float
);

template<ElementOrder O>
Eigen::Matrix<Float, NODES_IN_VOL_ELEM<O>, 1UL> shape_func_vol(
    const Float,
    const Float,
    const Float
);

template<ElementOrder O>
Eigen::Matrix<Float, DIM, NODES_IN_VOL_ELEM<O>> shape_func_vol_gradient(
    const Float,
    const Float,
    const Float
);

} // namespace numav
