// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/fem-helmholtz/shape-functions.hpp"

namespace numav {

template<>
Eigen::Matrix<Float, ENIS_COUNT<ElementOrder::LINEAR>, 1UL>
shape_func_sfc<ElementOrder::LINEAR>(
    const Float xi0,
    const Float xi1
) {
    const Float xi2 = 1_F - xi0 - xi1;
    return Eigen::Matrix<Float, ENIS_COUNT<ElementOrder::LINEAR>, 1UL>(
        xi0,
        xi1,
        xi2
    );
}

template<>
Eigen::Matrix<Float, ENIS_COUNT<ElementOrder::QUADRATIC>, 1UL>
shape_func_sfc<ElementOrder::QUADRATIC>(
    const Float xi0,
    const Float xi1
) {
    const Float xi2 = 1_F - xi0 - xi1;
    return Eigen::Matrix<Float, ENIS_COUNT<ElementOrder::QUADRATIC>, 1UL>(
        xi0*(2_F*xi0 - 1_F),
        xi1*(2_F*xi1 - 1_F),
        xi2*(2_F*xi2 - 1_F),
        4_F*xi0*xi1,
        4_F*xi0*xi2,
        4_F*xi1*xi2
    );
}

template<>
Eigen::Matrix<Float, ENIS_COUNT<ElementOrder::LINEAR>, 2UL>
shape_func_sfc_gradient<ElementOrder::LINEAR>(
    const Float xi0,
    const Float xi1
) {
    (void)xi0;
    (void)xi1;
    return Eigen::Matrix<Float, ENIS_COUNT<ElementOrder::LINEAR>, 2UL> {
        {+1_F, +0_F},
        {+0_F, +1_F},
        {-1_F, -1_F}
    };
}

template<>
Eigen::Matrix<Float, ENIS_COUNT<ElementOrder::QUADRATIC>, 2UL>
shape_func_sfc_gradient<ElementOrder::QUADRATIC>(
    const Float xi0,
    const Float xi1
) {
    const Float xi2 = 1_F - xi0 - xi1;
    return Eigen::Matrix<Float, ENIS_COUNT<ElementOrder::QUADRATIC>, 2UL> {
        { 4_F*xi0-1,        0_F             },
        { 0_F,              4_F*xi1-1       },
        { 1_F - 4_F*xi2,    1_F - 4_F*xi2   },
        { 4_F*xi1,          4_F*xi0         },
        { 4_F*(xi2 - xi0), -4_F*xi0         },
        {-4_F*xi1,          4_F*(xi2 - xi1) }
    };
}

template<>
Eigen::Matrix<Float, ENIV_COUNT<ElementOrder::LINEAR>, 1UL>
shape_func_vol<ElementOrder::LINEAR>(
    const Float xi0,
    const Float xi1,
    const Float xi2
) {
    const Float xi3 = 1_F - xi0 - xi1 - xi2;
    return Eigen::Matrix<Float, ENIV_COUNT<ElementOrder::LINEAR>, 1UL>(
        xi0,
        xi1,
        xi2,
        xi3
    );
}

template<>
Eigen::Matrix<Float, ENIV_COUNT<ElementOrder::QUADRATIC>, 1UL>
shape_func_vol<ElementOrder::QUADRATIC>(
    const Float xi0,
    const Float xi1,
    const Float xi2
) {
    const Float xi3 = 1_F - xi0 - xi1 - xi2;
    return Eigen::Matrix<Float, ENIV_COUNT<ElementOrder::QUADRATIC>, 1UL>(
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
Eigen::Matrix<Float, ENIV_COUNT<ElementOrder::LINEAR>, DIM>
shape_func_vol_gradient<ElementOrder::LINEAR>(
    const Float xi0,
    const Float xi1,
    const Float xi2
) {
    (void)xi0;
    (void)xi1;
    (void)xi2;
    return Eigen::Matrix<Float, ENIV_COUNT<ElementOrder::LINEAR>, DIM> {
        {+1_F, +0_F, +0_F},
        {+0_F, +1_F, +0_F},
        {+0_F, +0_F, +1_F},
        {-1_F, -1_F, -1_F}
    };
}

template<>
Eigen::Matrix<Float, ENIV_COUNT<ElementOrder::QUADRATIC>, DIM>
shape_func_vol_gradient<ElementOrder::QUADRATIC>(
    const Float xi0,
    const Float xi1,
    const Float xi2
) {
    const Float xi3 = 1_F - xi0 - xi1 - xi2;
    return Eigen::Matrix<Float, ENIV_COUNT<ElementOrder::QUADRATIC>, DIM> {
        { 4_F*xi0 - 1_F,    0_F,              0_F             },
        { 0_F,              4_F*xi1 - 1_F,    0_F             },
        { 0_F,              0_F,              4_F*xi2 - 1_F   },
        { 1_F - 4_F*xi3,    1_F - 4_F*xi3,    1_F - 4_F*xi3   },
        { 4_F*xi1,          4_F*xi0,          0_F             },
        { 4_F*xi2,          0_F,              4_F*xi0         },
        { 4_F*(xi3 - xi0), -4_F*xi0,         -4_F*xi0         },
        { 0_F,              4_F*xi2,          4_F*xi1         },
        {-4_F*xi1,          4_F*(xi3 - xi1), -4_F*xi1         },
        {-4_F*xi2,         -4_F*xi2,          4_F*(xi3 - xi2) }
    };
}

} // namespace numav
