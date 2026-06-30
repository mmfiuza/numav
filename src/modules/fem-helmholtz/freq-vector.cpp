// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/fem-helmholtz/impl.hpp"

#include "common/maths.hpp"

namespace numav {

template <ElementOrder O>
void SimulationFemHelmTet<O>::Impl::_define_freq_vector() {
    if(_freq_type_defined_by_user == _FreqTypeDefinedByUser::MAXIMUM) {
        _freq_min = 0_F;
    }
    if (_freq_type_defined_by_user != _FreqTypeDefinedByUser::STEPS) {
        switch(_freq_sampling_density) {
            case FrequencySamplingDensity::CONSTANT:
                _fi_to_freq = linspace(_freq_min, _freq_max, _fi_count);
                break;
            case FrequencySamplingDensity::QUADRATIC:
                _fi_to_freq = cubspace(_freq_min, _freq_max, _fi_count);
                break;
        }
    }
}

} // namespace numav

NUMAV_INSTANTIATE_SIM_AC_FEM_FREQ_D3
