// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/numav.hpp"

#include <array>

#include "Eigen/Eigen"

namespace numav {

// general get_stif_matrix_const_part declaration for all orders
template<ElementOrder O>
Eigen::Matrix<double,NODES_IN_VOL_ELEM<O>,NODES_IN_VOL_ELEM<O>>
get_stif_matrix_const_part(const std::array<std::array<double,3UL>,4UL>);

template<>
Eigen::Matrix<double,
    NODES_IN_VOL_ELEM<ElementOrder::O1>, NODES_IN_VOL_ELEM<ElementOrder::O1>
>
get_stif_matrix_const_part<ElementOrder::O1>(
    const std::array<std::array<double,3UL>,4UL> coords
) {
    const double x03 = coords[0UL][0UL] - coords[3UL][0UL];
    const double x13 = coords[1UL][0UL] - coords[3UL][0UL];
    const double x23 = coords[2UL][0UL] - coords[3UL][0UL];
    const double y03 = coords[0UL][1UL] - coords[3UL][1UL];
    const double y13 = coords[1UL][1UL] - coords[3UL][1UL];
    const double y23 = coords[2UL][1UL] - coords[3UL][1UL];
    const double z03 = coords[0UL][2UL] - coords[3UL][2UL];
    const double z13 = coords[1UL][2UL] - coords[3UL][2UL];
    const double z23 = coords[2UL][2UL] - coords[3UL][2UL];

    const double a = x23*z13 - x13*z23;
    const double b = y13*z23 - y23*z13;
    const double c = x13*y23 - x23*y13;
    const double d = x03*z23 - x23*z03;
    const double e = x23*y03 - x03*y23;
    const double f = y23*z03 - y03*z23;
    const double g = y03*z13 - y13*z03;
    const double h = x13*z03 - x03*z13;
    const double i = x03*y13 - x13*y03;
    const double j = -(c + i + e);
    const double k = -(f + b + g);
    const double l = -(d + a + h);

    Eigen::Matrix<double,
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
Eigen::Matrix<double,
    NODES_IN_VOL_ELEM<ElementOrder::O2>, NODES_IN_VOL_ELEM<ElementOrder::O2>
>
get_stif_matrix_const_part<ElementOrder::O2>(
    const std::array<std::array<double,3UL>,4UL> coords
) {
    // TODO
    error("Numav doesn't support second order analytic intgration for stifness"
        " matrices yet");

    (void) coords;

    Eigen::Matrix<double,
        NODES_IN_VOL_ELEM<ElementOrder::O2>,NODES_IN_VOL_ELEM<ElementOrder::O2>
    > btb;

    return btb;
}

template<ElementOrder O>
Eigen::Matrix<double,NODES_IN_VOL_ELEM<O>,NODES_IN_VOL_ELEM<O>> 
MASS_MATRIX_CONST_PART = [] {
    if constexpr (O == ElementOrder::O1) {
        return 
        Eigen::Matrix<double,NODES_IN_VOL_ELEM<O>,NODES_IN_VOL_ELEM<O>> {
            {+1.0/10.0, +1.0/20.0, +1.0/20.0, +1.0/20.0},
            {+1.0/20.0, +1.0/10.0, +1.0/20.0, +1.0/20.0},
            {+1.0/20.0, +1.0/20.0, +1.0/10.0, +1.0/20.0},
            {+1.0/20.0, +1.0/20.0, +1.0/20.0, +1.0/10.0}
        };
    }
    if constexpr (O == ElementOrder::O2) {
        return 
        Eigen::Matrix<double,NODES_IN_VOL_ELEM<O>,NODES_IN_VOL_ELEM<O>> {
            { +1.0/70.0, +1.0/420.0, +1.0/420.0, +1.0/420.0, -1.0/105.0,
                -1.0/105.0, -1.0/105.0,  -1.0/70.0,  -1.0/70.0,  -1.0/70.0},
            {+1.0/420.0,  +1.0/70.0, +1.0/420.0, +1.0/420.0, -1.0/105.0,
                 -1.0/70.0,  -1.0/70.0, -1.0/105.0, -1.0/105.0,  -1.0/70.0},
            {+1.0/420.0, +1.0/420.0,  +1.0/70.0, +1.0/420.0,  -1.0/70.0,
                -1.0/105.0,  -1.0/70.0, -1.0/105.0,  -1.0/70.0, -1.0/105.0},
            {+1.0/420.0, +1.0/420.0, +1.0/420.0,  +1.0/70.0,  -1.0/70.0,
                 -1.0/70.0, -1.0/105.0,  -1.0/70.0, -1.0/105.0, -1.0/105.0},
            {-1.0/105.0, -1.0/105.0,  -1.0/70.0,  -1.0/70.0, +8.0/105.0,
                +4.0/105.0, +4.0/105.0, +4.0/105.0, +4.0/105.0, +2.0/105.0},
            {-1.0/105.0,  -1.0/70.0, -1.0/105.0,  -1.0/70.0, +4.0/105.0,
                +8.0/105.0, +4.0/105.0, +4.0/105.0, +2.0/105.0, +4.0/105.0},
            {-1.0/105.0,  -1.0/70.0,  -1.0/70.0, -1.0/105.0, +4.0/105.0,
                +4.0/105.0, +8.0/105.0, +2.0/105.0, +4.0/105.0, +4.0/105.0},
            { -1.0/70.0, -1.0/105.0, -1.0/105.0,  -1.0/70.0, +4.0/105.0,
                +4.0/105.0, +2.0/105.0, +8.0/105.0, +4.0/105.0, +4.0/105.0},
            { -1.0/70.0, -1.0/105.0,  -1.0/70.0, -1.0/105.0, +4.0/105.0,
                +2.0/105.0, +4.0/105.0, +4.0/105.0, +8.0/105.0, +4.0/105.0},
            { -1.0/70.0,  -1.0/70.0, -1.0/105.0, -1.0/105.0, +2.0/105.0,
                +4.0/105.0, +4.0/105.0, +4.0/105.0, +4.0/105.0, +8.0/105.0}
        };
    }
}();

template<ElementOrder O>
Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,NODES_IN_SFC_ELEM<O>> 
DAMP_MATRIX_CONST_PART = [] {
    if constexpr (O == ElementOrder::O1) {
        return 
        Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,NODES_IN_SFC_ELEM<O>> {
            { +1.0/6.0, +1.0/12.0, +1.0/12.0},
            {+1.0/12.0,  +1.0/6.0, +1.0/12.0},
            {+1.0/12.0, +1.0/12.0,  +1.0/6.0}
        };
    }
    if constexpr (O == ElementOrder::O2) {
        return
        Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,NODES_IN_SFC_ELEM<O>> {
        { +1.0/30.0, -1.0/180.0, -1.0/180.0,      +0.0,      +0.0, -1.0/45.0},
        {-1.0/180.0,  +1.0/30.0, -1.0/180.0,      +0.0, -1.0/45.0,      +0.0},
        {-1.0/180.0, -1.0/180.0,  +1.0/30.0, -1.0/45.0,      +0.0,      +0.0},
        {      +0.0,       +0.0,  -1.0/45.0, +8.0/45.0, +4.0/45.0, +4.0/45.0},
        {      +0.0,  -1.0/45.0,       +0.0, +4.0/45.0, +8.0/45.0, +4.0/45.0},
        { -1.0/45.0,       +0.0,       +0.0, +4.0/45.0, +4.0/45.0, +8.0/45.0}
        };
    }
}();

template<ElementOrder O>
Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,1UL> 
FORC_VECTOR_CONST_PART = [] {
    if constexpr (O == ElementOrder::O1) {
        return Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,1UL> {
            1.0/3.0, 1.0/3.0, 1.0/3.0
        };
    }
    if constexpr (O == ElementOrder::O2) {
        return Eigen::Matrix<double,NODES_IN_SFC_ELEM<O>,1UL> {
            0.0, 0.0, 0.0, 1.0/3.0, 1.0/3.0, 1.0/3.0
        };
    }
}();

} // namespace numav