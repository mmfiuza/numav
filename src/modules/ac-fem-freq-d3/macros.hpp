// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

// explicit instantiation declarations
#define NUMAV_INSTANTIATE_SIM_AC_FEM_FREQ_D3 \
    namespace numav { \
        template class Simulation< \
            Phenomenon::ACOUSTIC, \
            NumericalMethod::FEM, \
            Domain::FREQUENCY, \
            Dimension::D3, \
            ElementShape::TETRAHEDRON, \
            ElementOrder::O1 \
        >; \
        template class Simulation< \
            Phenomenon::ACOUSTIC, \
            NumericalMethod::FEM, \
            Domain::FREQUENCY, \
            Dimension::D3, \
            ElementShape::TETRAHEDRON, \
            ElementOrder::O2 \
        >; \
    }
