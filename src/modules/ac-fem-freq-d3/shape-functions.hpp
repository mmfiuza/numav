// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/numav.hpp"
#include "modules/ac-fem-freq-d3/constants.hpp"
#include "Eigen/Eigen"

namespace numav {

template<ElementOrder O>
Eigen::Matrix<double, NODES_IN_SFC_ELEM<O>, 1> shape_func_sfc(
    const double&, const double&
);

template<ElementOrder O>
Eigen::Matrix<double, 2, NODES_IN_SFC_ELEM<O>> shape_func_sfc_gradient(
    const double&, const double&
);

template<ElementOrder O>
Eigen::Matrix<double, NODES_IN_VOL_ELEM<O>, 1> shape_func_vol(
    const double&, const double&, const double&
);

template<ElementOrder O>
Eigen::Matrix<double, DIM, NODES_IN_VOL_ELEM<O>> shape_func_vol_gradient(
    const double&, const double&, const double&
);

} // namespace numav
