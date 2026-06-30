// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/numav.hpp"
#include "H5Cpp.h"

namespace numav {

const H5std_string HDF5_SIGNATURE("numav_result_hdf5_0.2.0");

template<NumericalMethod N>
const H5std_string HDF5_NUMERICAL_METHOD = [] {
    if constexpr (N == NumericalMethod::FEM) {
        return "finite_element_method";
    }
}();

void write_string_attr(
    H5::H5Object& obj,
    const std::string& attr_name,
    const std::string& value
);

void write_uint64_attr(
    H5::H5Object& obj,
    const std::string& attr_name,
    const uint64_t value
);

H5::DataSet write_dataset_float64_1d(
    H5::Group& group,
    const std::string& name,
    const double* const data,
    const uint64_t size
);

H5::DataSet write_dataset_uint64_1d(
    H5::Group& group,
    const std::string& name,
    const uint64_t* const data,
    const uint64_t size
);

H5::DataSet write_dataset_float64_2d(
    H5::Group& group,
    const std::string& name,
    const double* const data,
    const uint64_t row_count,
    const uint64_t col_count
);

H5::DataSet write_dataset_uint64_2d(
    H5::Group& group,
    const std::string& name,
    const uint64_t* const data,
    const uint64_t row_count,
    const uint64_t col_count
);

H5::CompType get_hdf5_cmplx_type();

H5::DataSet write_dataset_cmplx128_2d(
    H5::Group& group,
    const std::string& name,
    const std::complex<double>* const data,
    const uint64_t row_count,
    const uint64_t col_count
);

H5::DataSet create_cmplx128_dataset(
    H5::Group& group,
    const std::string& name,
    const uint64_t row_count,
    const uint64_t col_count
);

} // namespace numav
