// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/numav.hpp"
#include "H5Cpp.h"

namespace numav {

const H5std_string HDF5_SIGNATURE("numav_result_hdf5_1.0.0");

template<Phenomenon P>
const H5std_string HDF5_PHENOMENON = [] {
    if constexpr (P == Phenomenon::ACOUSTIC) {
        return "acoustic";
    }
}();
template<NumericalMethod N>
const H5std_string HDF5_NUMERICAL_METHOD = [] {
    if constexpr (N == NumericalMethod::FEM) {
        return "finite_element_method";
    }
}();
template<Domain D>
const H5std_string HDF5_DOMAIN = [] {
    if constexpr (D == Domain::FREQUENCY) {
        return "frequency";
    }
    if constexpr (D == Domain::TIME) {
        return "time";
    }
}();
template<Dimension D>
const H5std_string HDF5_DIMENSION = [] {
    if constexpr (D == Dimension::D1) {
        return "1d";
    }
    if constexpr (D == Dimension::D2) {
        return "2d";
    }
    if constexpr (D == Dimension::D3) {
        return "3d";
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
    const size_t size
);

H5::DataSet write_dataset_uint64_1d(
    H5::Group& group,
    const std::string& name,
    const size_t* const data,
    const size_t size
);

H5::DataSet write_dataset_float64_2d(
    H5::Group& group,
    const std::string& name,
    const double* const data,
    const size_t row_count,
    const size_t col_count
);

H5::DataSet write_dataset_uint64_2d(
    H5::Group& group,
    const std::string& name,
    const size_t* const data,
    const size_t row_count,
    const size_t col_count
);

H5::CompType get_hdf5_cmplx_type();

H5::DataSet write_dataset_cmplx128_2d(
    H5::Group& group,
    const std::string& name,
    const std::complex<double>* const data,
    const size_t row_count,
    const size_t col_count
);

H5::DataSet create_cmplx128_dataset(
    H5::Group& group,
    const std::string& name,
    const size_t row_count,
    const size_t col_count
);

} // namespace numav
