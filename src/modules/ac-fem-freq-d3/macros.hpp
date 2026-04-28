// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

// explicit instantiation declarations
#define INSTANTIATE_SIMULATION_CLASS \
    template class Simulation< \
        Phenomenon::ACOUSTIC, \
        NumericalMethod::FEM, \
        Domain::FREQUENCY, \
        Dimension::D3, \
        ElementOrder::O1 \
    >; \
    template class Simulation< \
        Phenomenon::ACOUSTIC, \
        NumericalMethod::FEM, \
        Domain::FREQUENCY, \
        Dimension::D3, \
        ElementOrder::O2 \
    >;
