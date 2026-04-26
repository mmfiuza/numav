// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/shape-functions.hpp"

namespace numav {

template<>
Eigen::Matrix<Float, NODES_IN_SFC_ELEM<ElementOrder::O1>, 1UL>
shape_func_sfc<ElementOrder::O1>(
    const Float xi0,
    const Float xi1
) {
    const Float xi2 = 1.0 - xi0 - xi1;
    return Eigen::Matrix<Float, NODES_IN_SFC_ELEM<ElementOrder::O1>, 1UL>(
        xi0,
        xi1,
        xi2
    );
}

template<>
Eigen::Matrix<Float, NODES_IN_SFC_ELEM<ElementOrder::O2>, 1UL>
shape_func_sfc<ElementOrder::O2>(
    const Float xi0,
    const Float xi1
) {
    const Float xi2 = 1.0 - xi0 - xi1;
    return Eigen::Matrix<Float, NODES_IN_SFC_ELEM<ElementOrder::O2>, 1UL>(
        xi0*(2.0*xi0-1.0),
        xi1*(2.0*xi1-1.0),
        xi2*(2.0*xi2-1.0),
        4.0*xi0*xi1,
        4.0*xi0*xi2,
        4.0*xi1*xi2
    );
}

template<>
Eigen::Matrix<Float, 2UL, NODES_IN_SFC_ELEM<ElementOrder::O1>>
shape_func_sfc_gradient<ElementOrder::O1>(
    const Float xi0,
    const Float xi1
) {
    (void)xi0;
    (void)xi1;
    return Eigen::Matrix<Float, 2UL, NODES_IN_SFC_ELEM<ElementOrder::O1>> {
        {+1.0, +0.0, -1.0},
        {+0.0, +1.0, -1.0}
    };
}

template<>
Eigen::Matrix<Float, 2UL, NODES_IN_SFC_ELEM<ElementOrder::O2>>
shape_func_sfc_gradient<ElementOrder::O2>(
    const Float xi0,
    const Float xi1
) {
    const Float xi2 = 1.0 - xi0 - xi1;
    return Eigen::Matrix<Float, 2UL, NODES_IN_SFC_ELEM<ElementOrder::O2>> {
        {4.0*xi0-1, 0.0 , 1.0-4.0*xi2, 4.0*xi1, 4.0*(xi2-xi0), -4.0*xi1},
        {0.0, 4.0*xi1-1, 1.0-4.0*xi2, 4.0*xi0, -4.0*xi0, 4.0*(xi2-xi1)}
    };
}

template<>
Eigen::Matrix<Float, NODES_IN_VOL_ELEM<ElementOrder::O1>, 1UL>
shape_func_vol<ElementOrder::O1>(
    const Float xi0,
    const Float xi1,
    const Float xi2
) {
    const Float xi3 = 1.0 - xi0 - xi1 - xi2;
    return Eigen::Matrix<Float, NODES_IN_VOL_ELEM<ElementOrder::O1>, 1UL>(
        xi0,
        xi1,
        xi2,
        xi3
    );
}

template<>
Eigen::Matrix<Float, NODES_IN_VOL_ELEM<ElementOrder::O2>, 1UL>
shape_func_vol<ElementOrder::O2>(
    const Float xi0,
    const Float xi1,
    const Float xi2
) {
    const Float xi3 = 1.0 - xi0 - xi1 - xi2;
    return Eigen::Matrix<Float, NODES_IN_VOL_ELEM<ElementOrder::O2>, 1UL>(
        xi0*(2.0*xi0-1.0),
        xi1*(2.0*xi1-1.0),
        xi2*(2.0*xi2-1.0),
        xi3*(2.0*xi3-1.0),
        4.0*xi0*xi1,
        4.0*xi0*xi2,
        4.0*xi0*xi3,
        4.0*xi1*xi2,
        4.0*xi1*xi3,
        4.0*xi2*xi3
    );
}

template<>
Eigen::Matrix<Float, DIM, NODES_IN_VOL_ELEM<ElementOrder::O1>>
shape_func_vol_gradient<ElementOrder::O1>(
    const Float xi0,
    const Float xi1,
    const Float xi2
) {
    (void)xi0;
    (void)xi1;
    (void)xi2;
    return Eigen::Matrix<Float, DIM, NODES_IN_VOL_ELEM<ElementOrder::O1>> {
        {+1.0, +0.0, +0.0, -1.0},
        {+0.0, +1.0, +0.0, -1.0},
        {+0.0, +0.0, +1.0, -1.0}
    };
}

template<>
Eigen::Matrix<Float, DIM, NODES_IN_VOL_ELEM<ElementOrder::O2>>
shape_func_vol_gradient<ElementOrder::O2>(
    const Float xi0,
    const Float xi1,
    const Float xi2
) {
    const Float xi3 = 1.0 - xi0 - xi1 - xi2;
    return Eigen::Matrix<Float, DIM, NODES_IN_VOL_ELEM<ElementOrder::O2>> {
        {4.0*xi0-1.0, 0.0, 0.0, 1.0-4.0*xi3, 4.0*xi1,
            4.0*xi2, 4.0*(xi3-xi0), 0.0, -4.0*xi1, -4.0*xi2},
        {0.0, 4.0*xi1-1.0, 0.0, 1.0-4.0*xi3, 4.0*xi0,
            0.0, -4.0*xi0, 4.0*xi2, 4.0*(xi3-xi1), -4.0*xi2},
        {0.0, 0.0, 4.0*xi2-1.0, 1.0-4.0*xi3, 0.0,
            4.0*xi0, -4.0*xi0, 4.0*xi1, -4.0*xi1, 4.0*(xi3-xi2)}
    };
}

} // namespace numav
