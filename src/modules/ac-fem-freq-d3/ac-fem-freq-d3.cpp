// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/impl.hpp"

#define NUMAV_PARAMS_0()
#define NUMAV_PARAMS_1(t1) t1 v1
#define NUMAV_PARAMS_2(t1, t2) t1 v1, t2 v2
#define NUMAV_PARAMS_3(t1, t2, t3) t1 v1, t2 v2, t3 v3
#define NUMAV_PARAMS_4(t1,t2,t3,t4) t1 v1, t2 v2, t3 v3, t4 v4

#define NUMAV_ARGS_0()
#define NUMAV_ARGS_1(t1) v1
#define NUMAV_ARGS_2(t1, t2) v1, v2
#define NUMAV_ARGS_3(t1, t2, t3) v1, v2, v3
#define NUMAV_ARGS_4(t1,t2,t3,t4) v1, v2, v3, v4

#define NUMAV_COUNT(_0, _1, _2, _3, _4, N, ...) N
#define NUMAV_NARGS(...) NUMAV_COUNT(__VA_ARGS__ __VA_OPT__(,) 0,4,3,2,1,0)

#define NUMAV_CAT(a, b)  NUMAV_CAT_(a, b)
#define NUMAV_CAT_(a, b) a##b

#define NUMAV_PARAMS(...)\
    NUMAV_CAT(NUMAV_PARAMS_, NUMAV_NARGS(__VA_ARGS__))(__VA_ARGS__)
#define NUMAV_ARGS(...)\
    NUMAV_CAT(NUMAV_ARGS_, NUMAV_NARGS(__VA_ARGS__))(__VA_ARGS__)

#undef NUMAV_PUBLIC_METHOD
#define NUMAV_PUBLIC_METHOD(method_name, ...)\
    template <ElementOrder O>\
    void SimulationAcFemFreqD3<O>::method_name(NUMAV_PARAMS(__VA_ARGS__)) {\
        _pimpl->method_name(NUMAV_ARGS(__VA_ARGS__));\
    }

namespace numav {

NUMAV_SIM_AC_FEM_FREQ_D3_PUBLIC_METHODS

} // namespace numav

NUMAV_INSTANTIATE_SIM_AC_FEM_FREQ_D3
