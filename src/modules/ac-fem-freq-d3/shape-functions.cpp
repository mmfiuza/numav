// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/shape-functions.hpp"

namespace numav {

template<>
Eigen::Matrix<Float, NODES_IN_SFC_ELEM<ElementOrder::O1>, 1UL>
shape_func_sfc<ElementOrder::O1>(
    const Float xi0,
    const Float xi1
) {
    const Float xi2 = 1_F - xi0 - xi1;
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
    const Float xi2 = 1_F - xi0 - xi1;
    return Eigen::Matrix<Float, NODES_IN_SFC_ELEM<ElementOrder::O2>, 1UL>(
        xi0*(2_F*xi0 - 1_F),
        xi1*(2_F*xi1 - 1_F),
        xi2*(2_F*xi2 - 1_F),
        4_F*xi0*xi1,
        4_F*xi0*xi2,
        4_F*xi1*xi2
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
        {+1_F, +0_F, -1_F},
        {+0_F, +1_F, -1_F}
    };
}

template<>
Eigen::Matrix<Float, 2UL, NODES_IN_SFC_ELEM<ElementOrder::O2>>
shape_func_sfc_gradient<ElementOrder::O2>(
    const Float xi0,
    const Float xi1
) {
    const Float xi2 = 1_F - xi0 - xi1;
    return Eigen::Matrix<Float, 2UL, NODES_IN_SFC_ELEM<ElementOrder::O2>> {
        {
            4_F*xi0-1,
            0_F,
            1_F - 4_F*xi2,
            4_F*xi1,
            4_F*(xi2 - xi0),
            -4_F*xi1
        }, {
            0_F,
            4_F*xi1-1,
            1_F - 4_F*xi2,
            4_F*xi0,
            -4_F*xi0,
            4_F*(xi2 - xi1)
        }
    };
}

template<>
Eigen::Matrix<Float, NODES_IN_VOL_ELEM<ElementOrder::O1>, 1UL>
shape_func_vol<ElementOrder::O1>(
    const Float xi0,
    const Float xi1,
    const Float xi2
) {
    const Float xi3 = 1_F - xi0 - xi1 - xi2;
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
    const Float xi3 = 1_F - xi0 - xi1 - xi2;
    return Eigen::Matrix<Float, NODES_IN_VOL_ELEM<ElementOrder::O2>, 1UL>(
        xi0*(2_F*xi0 - 1_F),
        xi1*(2_F*xi1 - 1_F),
        xi2*(2_F*xi2 - 1_F),
        xi3*(2_F*xi3 - 1_F),
        4_F*xi0*xi1,
        4_F*xi0*xi2,
        4_F*xi0*xi3,
        4_F*xi1*xi2,
        4_F*xi1*xi3,
        4_F*xi2*xi3
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
        {+1_F, +0_F, +0_F, -1_F},
        {+0_F, +1_F, +0_F, -1_F},
        {+0_F, +0_F, +1_F, -1_F}
    };
}

template<>
Eigen::Matrix<Float, DIM, NODES_IN_VOL_ELEM<ElementOrder::O2>>
shape_func_vol_gradient<ElementOrder::O2>(
    const Float xi0,
    const Float xi1,
    const Float xi2
) {
    const Float xi3 = 1_F - xi0 - xi1 - xi2;
    return Eigen::Matrix<Float, DIM, NODES_IN_VOL_ELEM<ElementOrder::O2>> {
        {
            4_F*xi0 - 1_F,
            0_F,
            0_F,
            1_F - 4_F*xi3,
            4_F*xi1,
            4_F*xi2,
            4_F*(xi3 - xi0),
            0_F,
            -4_F*xi1,
            -4_F*xi2
        }, {
            0_F,
            4_F*xi1 - 1_F,
            0_F,
            1_F - 4_F*xi3,
            4_F*xi0,
            0_F,
            -4_F*xi0,
            4_F*xi2,
            4_F*(xi3 - xi1),
            -4_F*xi2
        }, {
            0_F,
            0_F,
            4_F*xi2 - 1_F,
            1_F - 4_F*xi3,
            0_F,
            4_F*xi0,
            -4_F*xi0,
            4_F*xi1,
            -4_F*xi1,
            4_F*(xi3 - xi2)
        }
    };
}

} // namespace numav
