// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/impl.hpp"

#include "common/maths.hpp"

namespace numav {

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_define_freq_vector() {
    if(_freq_type_defined_by_user == FreqTypeDefinedByUser::MAXIMUM) {
        _freq_min = 0;
    }
    if (_freq_type_defined_by_user != FreqTypeDefinedByUser::STEPS) {
        switch(_frequency_sampling_density) {
            case FrequencySamplingDensity::CONSTANT:
                _freq_steps = linspace(_freq_min, _freq_max, _freq_count);
                break;
            case FrequencySamplingDensity::QUADRATIC:
                _freq_steps = cubspace(_freq_min, _freq_max, _freq_count);
                break;
        }
    }
}

// explicit instantiation declarations
INSTANTIATE_SIMULATION_CLASS

} // namespace numav