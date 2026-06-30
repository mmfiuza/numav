// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

// explicit instantiation declarations
#define NUMAV_INSTANTIATE_SIM_AC_FEM_FREQ_D3 \
    namespace numav { \
        template class Simulation< \
            NumericalMethod::FEM, \
            Equation::HELMHOLTZ, \
            ElementShape::TETRAHEDRON, \
            ElementOrder::LINEAR \
        >; \
        template class Simulation< \
            NumericalMethod::FEM, \
            Equation::HELMHOLTZ, \
            ElementShape::TETRAHEDRON, \
            ElementOrder::QUADRATIC \
        >; \
    }
