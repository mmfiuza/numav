// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/hdf5.hpp"

namespace numav {

void write_string_attr(
    H5::H5Object& obj,
    const std::string& attr_name,
    const std::string& value
) {
    H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
    str_type.setCset(H5T_CSET_UTF8);
    H5::DataSpace scalar_space(H5S_SCALAR);
    H5::Attribute attr = obj.createAttribute(attr_name, str_type, scalar_space);
    attr.write(str_type, value);
}

void write_uint64_attr(
    H5::H5Object& obj,
    const std::string& attr_name,
    const uint64_t value
) {
    H5::DataSpace scalar_space(H5S_SCALAR);
    H5::Attribute attr = obj.createAttribute(
        attr_name, H5::PredType::STD_U64LE, scalar_space
    );
    attr.write(H5::PredType::NATIVE_UINT64, &value);
}

H5::DataSet write_dataset_float64_1d(
    H5::Group& group,
    const std::string& name,
    const double* const data,
    const size_t size
) {
    const hsize_t size_array[1UL] = { size };
    const H5::DataSpace data_space(1UL, size_array);
    H5::DataSet data_set = group.createDataSet(
        name, H5::PredType::IEEE_F64LE, data_space
    );
    data_set.write(data, H5::PredType::NATIVE_DOUBLE);
    return data_set;
}

H5::DataSet write_dataset_uint64_1d(
    H5::Group& group,
    const std::string& name,
    const size_t* const data,
    const size_t size
) {
    const hsize_t size_array[1UL] = { size };
    const H5::DataSpace data_space(1UL, size_array);
    H5::DataSet data_set = group.createDataSet(
        name, H5::PredType::STD_U64LE, data_space
    );
    data_set.write(data, H5::PredType::NATIVE_UINT64);
    return data_set;
}

H5::DataSet write_dataset_float64_2d(
    H5::Group& group,
    const std::string& name,
    const double* const data,
    const size_t row_count,
    const size_t col_count
) {
    const hsize_t size[2UL] = { row_count, col_count };
    const H5::DataSpace data_space(2UL, size);
    H5::DataSet data_set = group.createDataSet(
        name, H5::PredType::IEEE_F64LE, data_space
    );
    data_set.write(data, H5::PredType::NATIVE_DOUBLE);
    return data_set;
}

H5::DataSet write_dataset_uint64_2d(
    H5::Group& group,
    const std::string& name,
    const size_t* const data,
    const size_t row_count,
    const size_t col_count
) {
    const hsize_t size[2UL] = { row_count, col_count };
    const H5::DataSpace data_space(2UL, size);
    H5::DataSet data_set = group.createDataSet(
        name, H5::PredType::STD_U64LE, data_space
    );
    data_set.write(data, H5::PredType::NATIVE_UINT64);
    return data_set;
}

H5::CompType get_hdf5_cmplx_type() {
    H5::CompType data_type(sizeof(std::complex<double>));
    data_type.insertMember("r", 0UL, H5::PredType::NATIVE_DOUBLE);
    data_type.insertMember("i", sizeof(double), H5::PredType::NATIVE_DOUBLE);
    return data_type;
}

H5::DataSet write_dataset_cmplx128_2d(
    H5::Group& group,
    const std::string& name,
    const std::complex<double>* const data,
    const size_t row_count,
    const size_t col_count
) {
    H5::CompType cmplx_type = get_hdf5_cmplx_type();
    const hsize_t size[2UL] = { row_count, col_count };
    const H5::DataSpace data_space(2UL, size);
    H5::DataSet data_set = group.createDataSet(name, cmplx_type, data_space);
    data_set.write(data, cmplx_type);
    return data_set;
}

H5::DataSet create_cmplx128_dataset(
    H5::Group& group,
    const std::string& name,
    const size_t row_count,
    const size_t col_count
) {
    H5::CompType cmplx_type = get_hdf5_cmplx_type();
    const hsize_t dims[2UL] = { row_count, col_count };

    // chunking
    const hsize_t chunk[2UL] = { 1UL, col_count }; // one freq row per chunk
    H5::DSetCreatPropList props;
    props.setChunk(2UL, chunk);
    props.setDeflate(6); // compression level

    const H5::DataSpace data_space(2UL, dims);
    return group.createDataSet(name, cmplx_type, data_space, props);
}

} // namespace numav
