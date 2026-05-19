// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "common/maths.hpp"
#include "modules/ac-fem-freq-d3/impl.hpp"

namespace numav {

template<ElementOrder O> constexpr size_t NGP_FORC = [] {
    if constexpr (O == ElementOrder::O1) { return 1UL; }
    if constexpr (O == ElementOrder::O2) { return 3UL; }
}();

template<ElementOrder O> constexpr size_t NGP_DAMP = [] {
    if constexpr (O == ElementOrder::O1) { return 1UL; }
    if constexpr (O == ElementOrder::O2) { return 3UL; }
}();

template<ElementOrder O> constexpr size_t NGP_STIF = [] {
    if constexpr (O == ElementOrder::O1) { return 1UL; }
    if constexpr (O == ElementOrder::O2) { return 4UL; }
}();

template<ElementOrder O> constexpr size_t NGP_MASS = [] {
    if constexpr (O == ElementOrder::O1) { return 4UL;  }
    if constexpr (O == ElementOrder::O2) { return 15UL; }
}();

template<size_t N>
constexpr std::array<std::array<Float,2UL>,N> GAUSS_POINTS_SFC = [] {
    if constexpr (N == 1UL) {
        constexpr Float a = 1_F / 3_F;
        return std::array<std::array<Float,2UL>,N>{{ {a, a} }};
    }
    if constexpr (N == 3UL) {
        constexpr Float a = 1_F / 6_F;
        constexpr Float b = 2_F / 3_F;
        return std::array<std::array<Float,2UL>,N>{{
            {a,a}, {b,a}, {a,b}
        }};
    }
}();

template<size_t N>
constexpr std::array<std::array<Float,DIM>,N> GAUSS_POINTS_VOL = [] {
    if constexpr (N == 1UL) {
        constexpr Float a = 1_F / 4_F;
        return std::array<std::array<Float,DIM>,N>{{ {a, a, a} }};
    }
    if constexpr (N == 4UL) {
        constexpr Float a = (5_F -     SQRT(5_F)) / 20_F;
        constexpr Float b = (5_F + 3_F*SQRT(5_F)) / 20_F;
        return std::array<std::array<Float,DIM>,N>{{
            {a,a,a}, {b,a,a}, {a,b,a}, {a,a,b}
        }};
    }
    if constexpr (N == 5UL) {
        constexpr Float a = 1_F / 4_F;
        constexpr Float b = 1_F / 6_F;
        constexpr Float c = 1_F / 2_F;
        return std::array<std::array<Float,DIM>,N>{{
            {a,a,a},
            {b,b,b}, {c,b,b}, {b,c,b}, {b,b,c}
        }};
    }
    if constexpr (N == 11UL) {
        constexpr Float a =  1_F / 4_F;
        constexpr Float b =  1_F / 14_F;
        constexpr Float c = 11_F / 14_F;
        constexpr Float d = 0.399403576166799_F; // TODO: discover the fraction
        constexpr Float e = 0.100596423833201_F; // TODO: discover the fraction
        return std::array<std::array<Float,DIM>,N>{{
            {a,a,a},
            {b,b,b}, {c,b,b}, {b,c,b}, {b,b,c},
            {d,d,e}, {d,e,d}, {e,d,d}, {d,e,e}, {e,d,e}, {e,e,d}
        }};
    }
    if constexpr (N == 15UL) {
        constexpr Float a = 1_F / 4_F;
        constexpr Float b = ( 7_F +     SQRT(15_F)) / 34_F;
        constexpr Float c = ( 7_F -     SQRT(15_F)) / 34_F;
        constexpr Float d = (13_F - 3_F*SQRT(15_F)) / 34_F;
        constexpr Float e = (13_F + 3_F*SQRT(15_F)) / 34_F;
        constexpr Float f = ( 5_F -     SQRT(15_F)) / 20_F;
        constexpr Float g = ( 5_F +     SQRT(15_F)) / 20_F;
        return std::array<std::array<Float,DIM>,N>{{
            {a,a,a},
            {b,b,b}, {b,b,d}, {b,d,b}, {d,b,b}, 
            {c,c,c}, {c,c,e}, {c,e,c}, {e,c,c},
            {f,f,g}, {f,g,f}, {g,f,f}, {f,g,g}, {g,f,g}, {g,g,f}
        }};
    }
}();

template<size_t N>
constexpr std::array<Float,N> GAUSS_WEIGHTS_SFC = [] {
    if constexpr (N == 1UL) {
        constexpr Float a = 1_F / 2_F;
        return std::array<Float,N>({a});
    }
    if constexpr (N == 3UL) {
        constexpr Float a = 1_F / 6_F;
        return std::array<Float,N>({a,a,a});
    }
}();

template<size_t N>
constexpr std::array<Float,N> GAUSS_WEIGHTS_VOL = [] {
    if constexpr (N == 1UL) {
        constexpr Float a = 1_F / 6_F;
        return std::array<Float,N>({a});
    }
    if constexpr (N == 4UL) {
        constexpr Float a = 1_F / 24_F;
        return std::array<Float,N>({a,a,a,a});
    }
    if constexpr (N == 5UL) {
        constexpr Float a = -2_F / 15_F;
        constexpr Float b = 3_F / 40_F;
        return std::array<Float,N>({a,b,b,b,b});
    }
    if constexpr (N == 11UL) {
        constexpr Float a = -74_F / 5625_F;
        constexpr Float b = 343_F / 45000_F;
        constexpr Float c = 28_F / 1125_F;
        return std::array<Float,N>({a,b,b,b,b,c,c,c,c,c,c});
    }
    if constexpr (N == 15UL) {
        constexpr Float a = 8_F / 405_F;
        constexpr Float b = (2665_F - 14_F*SQRT(15_F)) / 226800_F;
        constexpr Float c = (2665_F + 14_F*SQRT(15_F)) / 226800_F;
        constexpr Float d = 5_F / 567_F;
        return std::array<Float,N>({a,b,b,b,b,c,c,c,c,d,d,d,d,d,d});
    }
}();

} // namespace numav
