// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/fem-helmholtz/impl.hpp"

namespace numav {

template <ElementOrder O>
void SimulationFemHelmTet<O>::Impl::_post_process()
{
    // TODO   
}

} // namespace numav

NUMAV_INSTANTIATE_SIM_AC_FEM_FREQ_D3
