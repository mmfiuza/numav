// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/numav.hpp"

#include <array>

#include "Eigen/Eigen"

namespace numav {

// general get_stif_matrix_const_part declaration for all orders
template<ElementOrder O>
Eigen::Matrix<Float,NODES_IN_VOL_ELEM<O>,NODES_IN_VOL_ELEM<O>>
get_stif_matrix_const_part(const std::array<std::array<Float,3UL>,4UL>);

template<>
Eigen::Matrix<Float,
    NODES_IN_VOL_ELEM<ElementOrder::O1>, NODES_IN_VOL_ELEM<ElementOrder::O1>
>
get_stif_matrix_const_part<ElementOrder::O1>(
    const std::array<std::array<Float,3UL>,4UL> coords
) {
    const Float x03 = coords[0UL][0UL] - coords[3UL][0UL];
    const Float x13 = coords[1UL][0UL] - coords[3UL][0UL];
    const Float x23 = coords[2UL][0UL] - coords[3UL][0UL];
    const Float y03 = coords[0UL][1UL] - coords[3UL][1UL];
    const Float y13 = coords[1UL][1UL] - coords[3UL][1UL];
    const Float y23 = coords[2UL][1UL] - coords[3UL][1UL];
    const Float z03 = coords[0UL][2UL] - coords[3UL][2UL];
    const Float z13 = coords[1UL][2UL] - coords[3UL][2UL];
    const Float z23 = coords[2UL][2UL] - coords[3UL][2UL];

    const Float a = x23*z13 - x13*z23;
    const Float b = y13*z23 - y23*z13;
    const Float c = x13*y23 - x23*y13;
    const Float d = x03*z23 - x23*z03;
    const Float e = x23*y03 - x03*y23;
    const Float f = y23*z03 - y03*z23;
    const Float g = y03*z13 - y13*z03;
    const Float h = x13*z03 - x03*z13;
    const Float i = x03*y13 - x13*y03;
    const Float j = -(c + i + e);
    const Float k = -(f + b + g);
    const Float l = -(d + a + h);

    Eigen::Matrix<Float,
        NODES_IN_VOL_ELEM<ElementOrder::O1>,NODES_IN_VOL_ELEM<ElementOrder::O1>
    > btb;

    btb(0UL,0UL) = a*a + b*b + c*c;
    btb(0UL,1UL) = a*d + b*f + c*e;
    btb(0UL,2UL) = a*h + b*g + c*i;
    btb(0UL,3UL) = a*l + b*k + c*j;
    btb(1UL,1UL) = d*d + e*e + f*f;
    btb(1UL,2UL) = d*h + e*i + f*g;
    btb(1UL,3UL) = d*l + e*j + f*k;
    btb(2UL,2UL) = g*g + h*h + i*i;
    btb(2UL,3UL) = g*k + h*l + i*j;
    btb(3UL,3UL) = j*j + k*k + l*l;

    return btb;
}

template<>
Eigen::Matrix<Float,
    NODES_IN_VOL_ELEM<ElementOrder::O2>, NODES_IN_VOL_ELEM<ElementOrder::O2>
>
get_stif_matrix_const_part<ElementOrder::O2>(
    const std::array<std::array<Float,3UL>,4UL> coords
) {
    // TODO
    error("Numav doesn't support second order analytic intgration for stifness"
        " matrices yet");

    (void) coords;

    Eigen::Matrix<Float,
        NODES_IN_VOL_ELEM<ElementOrder::O2>,NODES_IN_VOL_ELEM<ElementOrder::O2>
    > btb;

    return btb;
}

template<ElementOrder O>
Eigen::Matrix<Float,NODES_IN_VOL_ELEM<O>,NODES_IN_VOL_ELEM<O>> 
MASS_MATRIX_CONST_PART = [] {
    if constexpr (O == ElementOrder::O1) {
        return 
        Eigen::Matrix<Float,NODES_IN_VOL_ELEM<O>,NODES_IN_VOL_ELEM<O>> {
            {+1_F/10_F, +1_F/20_F, +1_F/20_F, +1_F/20_F},
            {+1_F/20_F, +1_F/10_F, +1_F/20_F, +1_F/20_F},
            {+1_F/20_F, +1_F/20_F, +1_F/10_F, +1_F/20_F},
            {+1_F/20_F, +1_F/20_F, +1_F/20_F, +1_F/10_F}
        };
    }
    if constexpr (O == ElementOrder::O2) {
        return 
        Eigen::Matrix<Float,NODES_IN_VOL_ELEM<O>,NODES_IN_VOL_ELEM<O>> {
            { +1_F/70_F, +1_F/420_F, +1_F/420_F, +1_F/420_F, -1_F/105_F,
                -1_F/105_F, -1_F/105_F,  -1_F/70_F,  -1_F/70_F,  -1_F/70_F},
            {+1_F/420_F,  +1_F/70_F, +1_F/420_F, +1_F/420_F, -1_F/105_F,
                 -1_F/70_F,  -1_F/70_F, -1_F/105_F, -1_F/105_F,  -1_F/70_F},
            {+1_F/420_F, +1_F/420_F,  +1_F/70_F, +1_F/420_F,  -1_F/70_F,
                -1_F/105_F,  -1_F/70_F, -1_F/105_F,  -1_F/70_F, -1_F/105_F},
            {+1_F/420_F, +1_F/420_F, +1_F/420_F,  +1_F/70_F,  -1_F/70_F,
                 -1_F/70_F, -1_F/105_F,  -1_F/70_F, -1_F/105_F, -1_F/105_F},
            {-1_F/105_F, -1_F/105_F,  -1_F/70_F,  -1_F/70_F, +8_F/105_F,
                +4_F/105_F, +4_F/105_F, +4_F/105_F, +4_F/105_F, +2_F/105_F},
            {-1_F/105_F,  -1_F/70_F, -1_F/105_F,  -1_F/70_F, +4_F/105_F,
                +8_F/105_F, +4_F/105_F, +4_F/105_F, +2_F/105_F, +4_F/105_F},
            {-1_F/105_F,  -1_F/70_F,  -1_F/70_F, -1_F/105_F, +4_F/105_F,
                +4_F/105_F, +8_F/105_F, +2_F/105_F, +4_F/105_F, +4_F/105_F},
            { -1_F/70_F, -1_F/105_F, -1_F/105_F,  -1_F/70_F, +4_F/105_F,
                +4_F/105_F, +2_F/105_F, +8_F/105_F, +4_F/105_F, +4_F/105_F},
            { -1_F/70_F, -1_F/105_F,  -1_F/70_F, -1_F/105_F, +4_F/105_F,
                +2_F/105_F, +4_F/105_F, +4_F/105_F, +8_F/105_F, +4_F/105_F},
            { -1_F/70_F,  -1_F/70_F, -1_F/105_F, -1_F/105_F, +2_F/105_F,
                +4_F/105_F, +4_F/105_F, +4_F/105_F, +4_F/105_F, +8_F/105_F}
        };
    }
}();

template<ElementOrder O>
Eigen::Matrix<Float,NODES_IN_SFC_ELEM<O>,NODES_IN_SFC_ELEM<O>> 
DAMP_MATRIX_CONST_PART = [] {
    if constexpr (O == ElementOrder::O1) {
        return 
        Eigen::Matrix<Float,NODES_IN_SFC_ELEM<O>,NODES_IN_SFC_ELEM<O>> {
            { +1_F/6_F, +1_F/12_F, +1_F/12_F},
            {+1_F/12_F,  +1_F/6_F, +1_F/12_F},
            {+1_F/12_F, +1_F/12_F,  +1_F/6_F}
        };
    }
    if constexpr (O == ElementOrder::O2) {
        return
        Eigen::Matrix<Float,NODES_IN_SFC_ELEM<O>,NODES_IN_SFC_ELEM<O>> {
        { +1_F/30_F, -1_F/180_F, -1_F/180_F,      +0_F,      +0_F, -1_F/45_F},
        {-1_F/180_F,  +1_F/30_F, -1_F/180_F,      +0_F, -1_F/45_F,      +0_F},
        {-1_F/180_F, -1_F/180_F,  +1_F/30_F, -1_F/45_F,      +0_F,      +0_F},
        {      +0_F,       +0_F,  -1_F/45_F, +8_F/45_F, +4_F/45_F, +4_F/45_F},
        {      +0_F,  -1_F/45_F,       +0_F, +4_F/45_F, +8_F/45_F, +4_F/45_F},
        { -1_F/45_F,       +0_F,       +0_F, +4_F/45_F, +4_F/45_F, +8_F/45_F}
        };
    }
}();

template<ElementOrder O>
Eigen::Matrix<Float,NODES_IN_SFC_ELEM<O>,1UL> 
FORC_VECTOR_CONST_PART = [] {
    if constexpr (O == ElementOrder::O1) {
        return Eigen::Matrix<Float,NODES_IN_SFC_ELEM<O>,1UL> {
            1_F/3_F, 1_F/3_F, 1_F/3_F
        };
    }
    if constexpr (O == ElementOrder::O2) {
        return Eigen::Matrix<Float,NODES_IN_SFC_ELEM<O>,1UL> {
            0_F, 0_F, 0_F, 1_F/3_F, 1_F/3_F, 1_F/3_F
        };
    }
}();

} // namespace numav