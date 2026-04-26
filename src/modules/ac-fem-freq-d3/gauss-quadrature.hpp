// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

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
constexpr std::array<std::array<double,2UL>,N> GAUSS_POINTS_SFC = [] {
    if constexpr (N == 1UL) {
        constexpr double a = 1.0 / 3.0;
        return std::array<std::array<double,2UL>,N>({ {a, a} });
    }
    if constexpr (N == 3UL) {
        constexpr double a = 1.0 / 6.0;
        constexpr double b = 2.0 / 3.0;
        return std::array<std::array<double,2UL>,N>{{
            {a,a}, {b,a}, {a,b}
        }};
    }
}();

template<size_t N>
constexpr std::array<std::array<double,DIM>,N> GAUSS_POINTS_VOL = [] {
    if constexpr (N == 1UL) {
        constexpr double a = 1.0 / 4.0;
        return std::array<std::array<double,DIM>,N>({ {a, a, a} });
    }
    if constexpr (N == 4UL) {
        constexpr double a = (5.0 -     std::sqrt(5.0)) / 20.0;
        constexpr double b = (5.0 + 3.0*std::sqrt(5.0)) / 20.0;
        return std::array<std::array<double,DIM>,N>{{
            {a,a,a}, {b,a,a}, {a,b,a}, {a,a,b}
        }};
    }
    if constexpr (N == 5UL) {
        constexpr double a = 1.0 / 4.0;
        constexpr double b = 1.0 / 6.0;
        constexpr double c = 1.0 / 2.0;
        return std::array<std::array<double,DIM>,N>{{
            {a,a,a},
            {b,b,b}, {c,b,b}, {b,c,b}, {b,b,c}
        }};
    }
    if constexpr (N == 11UL) {
        constexpr double a =  1.0 / 4.0;
        constexpr double b =  1.0 / 14.0;
        constexpr double c = 11.0 / 14.0;
        constexpr double d = 0.399403576166799; // TODO: discover the fraction
        constexpr double e = 0.100596423833201; // TODO: discover the fraction
        return std::array<std::array<double,DIM>,N>{{
            {a,a,a},
            {b,b,b}, {c,b,b}, {b,c,b}, {b,b,c},
            {d,d,e}, {d,e,d}, {e,d,d}, {d,e,e}, {e,d,e}, {e,e,d}
        }};
    }
    if constexpr (N == 15UL) {
        constexpr double a = 1.0 / 4.0;
        constexpr double b = ( 7.0 +     std::sqrt(15.0)) / 34.0;
        constexpr double c = ( 7.0 -     std::sqrt(15.0)) / 34.0;
        constexpr double d = (13.0 - 3.0*std::sqrt(15.0)) / 34.0;
        constexpr double e = (13.0 + 3.0*std::sqrt(15.0)) / 34.0;
        constexpr double f = ( 5.0 -     std::sqrt(15.0)) / 20.0;
        constexpr double g = ( 5.0 +     std::sqrt(15.0)) / 20.0;
        return std::array<std::array<double,DIM>,N>{{
            {a,a,a},
            {b,b,b}, {b,b,d}, {b,d,b}, {d,b,b}, 
            {c,c,c}, {c,c,e}, {c,e,c}, {e,c,c},
            {f,f,g}, {f,g,f}, {g,f,f}, {f,g,g}, {g,f,g}, {g,g,f}
        }};
    }
}();

template<size_t N>
constexpr std::array<double,N> GAUSS_WEIGHTS_SFC = [] {
    if constexpr (N == 1UL) {
        constexpr double a = 1.0 / 2.0;
        return std::array<double,N>({a});
    }
    if constexpr (N == 3UL) {
        constexpr double a = 1.0 / 6.0;
        return std::array<double,N>({a,a,a});
    }
}();

template<size_t N>
constexpr std::array<double,N> GAUSS_WEIGHTS_VOL = [] {
    if constexpr (N == 1UL) {
        constexpr double a = 1.0 / 6.0;
        return std::array<double,N>({a});
    }
    if constexpr (N == 4UL) {
        constexpr double a = 1.0 / 24.0;
        return std::array<double,N>({a,a,a,a});
    }
    if constexpr (N == 5UL) {
        constexpr double a = -2.0 / 15.0;
        constexpr double b = 3.0 / 40.0;
        return std::array<double,N>({a,b,b,b,b});
    }
    if constexpr (N == 11UL) {
        constexpr double a = -74.0 / 5625.0;
        constexpr double b = 343.0 / 45000.0;
        constexpr double c = 28.0 / 1125.0;
        return std::array<double,N>({a,b,b,b,b,c,c,c,c,c,c});
    }
    if constexpr (N == 15UL) {
        constexpr double a = 8.0 / 405.0;
        constexpr double b = (2665.0 - 14.0*std::sqrt(15.0)) / 226800.0;
        constexpr double c = (2665.0 + 14.0*std::sqrt(15.0)) / 226800.0;
        constexpr double d = 5.0 / 567.0;
        return std::array<double,N>({a,b,b,b,b,c,c,c,c,d,d,d,d,d,d});
    }
}();

} // namespace numav
