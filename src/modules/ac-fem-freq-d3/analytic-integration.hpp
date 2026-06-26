// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/numav.hpp"
#include "common/maths.hpp"
#include "modules/ac-fem-freq-d3/shape-functions.hpp"

#include <array>

#include "Eigen/Eigen"

namespace numav {

// general get_stif_matrix_const_part declaration for all orders
template<ElementOrder O>
Eigen::Matrix<Float,ENIV_COUNT<O>,ENIV_COUNT<O>> get_stif_matrix_const_part(
    const Eigen::Matrix<Float,DIM,ENIV_COUNT<O>>
);

template<>
Eigen::Matrix<Float,ENIV_COUNT<ElementOrder::O1>,ENIV_COUNT<ElementOrder::O1>>
get_stif_matrix_const_part<ElementOrder::O1>(
    const Eigen::Matrix<Float,DIM,ENIV_COUNT<ElementOrder::O1>> coords
) {
    const Float x03 = coords(0UL, 0UL) - coords(0UL, 3UL);
    const Float x13 = coords(0UL, 1UL) - coords(0UL, 3UL);
    const Float x23 = coords(0UL, 2UL) - coords(0UL, 3UL);
    const Float y03 = coords(1UL, 0UL) - coords(1UL, 3UL);
    const Float y13 = coords(1UL, 1UL) - coords(1UL, 3UL);
    const Float y23 = coords(1UL, 2UL) - coords(1UL, 3UL);
    const Float z03 = coords(2UL, 0UL) - coords(2UL, 3UL);
    const Float z13 = coords(2UL, 1UL) - coords(2UL, 3UL);
    const Float z23 = coords(2UL, 2UL) - coords(2UL, 3UL);

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

    return Eigen::Matrix<
        Float,ENIV_COUNT<ElementOrder::O1>,ENIV_COUNT<ElementOrder::O1>
    > {
        {a*a + b*b + c*c, a*d + b*f + c*e, a*h + b*g + c*i, a*l + b*k + c*j},
        {              0, d*d + e*e + f*f, d*h + e*i + f*g, d*l + e*j + f*k},
        {              0,               0, g*g + h*h + i*i, g*k + h*l + i*j},
        {              0,               0,               0, j*j + k*k + l*l}
    } / 36_F;
}

Eigen::Matrix<Float,3UL,3UL> adjugate(const Eigen::Matrix<Float,3UL,3UL> mat) {
    const Float& a = mat(0,0);
    const Float& b = mat(0,1);
    const Float& c = mat(0,2);
    const Float& d = mat(1,0);
    const Float& e = mat(1,1);
    const Float& f = mat(1,2);
    const Float& g = mat(2,0);
    const Float& h = mat(2,1);
    const Float& i = mat(2,2);
    const Float cf11 = +(e*i - f*h);
    const Float cf12 = -(d*i - f*g);
    const Float cf13 = +(d*h - e*g);
    const Float cf21 = -(b*i - c*h);
    const Float cf22 = +(a*i - c*g);
    const Float cf23 = -(a*h - b*g);
    const Float cf31 = +(b*f - c*e);
    const Float cf32 = -(a*f - c*d);
    const Float cf33 = +(a*e - b*d);
    return Eigen::Matrix<Float,3UL,3UL> {
        {cf11, cf21, cf31},
        {cf12, cf22, cf32},
        {cf13, cf23, cf33}
    };
}

template<>
Eigen::Matrix<Float,ENIV_COUNT<ElementOrder::O2>,ENIV_COUNT<ElementOrder::O2>>
get_stif_matrix_const_part<ElementOrder::O2>(
    const Eigen::Matrix<Float,DIM,ENIV_COUNT<ElementOrder::O2>> coords
) {
    Eigen::Matrix<Float,
        ENIV_COUNT<ElementOrder::O2>, ENIV_COUNT<ElementOrder::O2>
    > matrix;
    matrix.setZero();

    const Float a = (5_F -     SQRT(5_F)) / 20_F;
    const Float b = (5_F + 3_F*SQRT(5_F)) / 20_F;
    const std::array<std::array<Float,3UL>,4UL> points = {{
        {a,a,a}, {b,a,a}, {a,b,a}, {a,a,b}
    }};
    for (size_t gpi = 0UL; gpi != 4UL; ++gpi)
    {
        const Eigen::Matrix<Float, ENIV_COUNT<ElementOrder::O2>, DIM> 
            nabla_n = shape_func_vol_gradient<ElementOrder::O2>(
                points[gpi][0UL], points[gpi][1UL], points[gpi][2UL]
            );
        
        const Eigen::Matrix<Float, ENIV_COUNT<ElementOrder::O2>, DIM> b =
            nabla_n * adjugate(coords * nabla_n);
        
        matrix += b * b.transpose();
    }
    matrix *= (1_F / 144_F);
    return matrix;
}

template<ElementOrder O>
Eigen::Matrix<Float,ENIV_COUNT<O>,ENIV_COUNT<O>> MASS_MATRIX_CONST_PART = [] {
    if constexpr (O == ElementOrder::O1) {
        return 
        Eigen::Matrix<Float,ENIV_COUNT<O>,ENIV_COUNT<O>> {
            {+1_F/10_F, +1_F/20_F, +1_F/20_F, +1_F/20_F},
            {+1_F/20_F, +1_F/10_F, +1_F/20_F, +1_F/20_F},
            {+1_F/20_F, +1_F/20_F, +1_F/10_F, +1_F/20_F},
            {+1_F/20_F, +1_F/20_F, +1_F/20_F, +1_F/10_F}
        };
    }
    if constexpr (O == ElementOrder::O2) {
        return 
        Eigen::Matrix<Float,ENIV_COUNT<O>,ENIV_COUNT<O>> {
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
Eigen::Matrix<Float,ENIS_COUNT<O>,ENIS_COUNT<O>> DAMP_MATRIX_CONST_PART = [] {
    if constexpr (O == ElementOrder::O1) {
        return 
        Eigen::Matrix<Float,ENIS_COUNT<O>,ENIS_COUNT<O>> {
            { +1_F/6_F, +1_F/12_F, +1_F/12_F},
            {+1_F/12_F,  +1_F/6_F, +1_F/12_F},
            {+1_F/12_F, +1_F/12_F,  +1_F/6_F}
        };
    }
    if constexpr (O == ElementOrder::O2) {
        return
        Eigen::Matrix<Float,ENIS_COUNT<O>,ENIS_COUNT<O>> {
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
Eigen::Matrix<Float,ENIS_COUNT<O>,1UL> FORC_VECTOR_CONST_PART = [] {
    if constexpr (O == ElementOrder::O1) {
        return Eigen::Matrix<Float,ENIS_COUNT<O>,1UL> {
            1_F/3_F, 1_F/3_F, 1_F/3_F
        };
    }
    if constexpr (O == ElementOrder::O2) {
        return Eigen::Matrix<Float,ENIS_COUNT<O>,1UL> {
            0_F, 0_F, 0_F, 1_F/3_F, 1_F/3_F, 1_F/3_F
        };
    }
}();

} // namespace numav{