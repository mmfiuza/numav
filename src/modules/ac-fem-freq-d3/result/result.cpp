// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav/numav.hpp"

namespace numav {

ResultAcFemFreqD3::Result() = default;
ResultAcFemFreqD3::Result(const size_t& node_count, const size_t& freq_count) {
    _data = Eigen::Matrix<_cmplx_t, Eigen::Dynamic, Eigen::Dynamic>(
        node_count, freq_count
    );
};
ResultAcFemFreqD3::~Result() = default;

} // namespace numav
