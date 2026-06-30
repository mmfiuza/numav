// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include <memory>
#include <array>
#include <vector>

#include "numav/aliases.hpp"

namespace numav {

    enum class NumericalMethod : uint64_t {
        FEM
    };

    // declare the general Simulation class
    template<NumericalMethod N, auto... ARGS> class Simulation {};

} // namespace numav

// include modules
#include "numav/modules/fem-helmholtz.hpp"
