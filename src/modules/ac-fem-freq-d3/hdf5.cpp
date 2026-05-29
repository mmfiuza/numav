// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/exception.hpp"
#include "common/utils.hpp"
#include "common/hdf5.hpp"
#include "modules/ac-fem-freq-d3/impl.hpp"

#include "H5DSpublic.h"

namespace numav {

template<ElementOrder O>
const H5std_string HDF5_ORDER = [] {
    if constexpr (O == ElementOrder::O1) {
        return "linear";
    }
    if constexpr (O == ElementOrder::O2) {
        return "quadratic";
    }
}();

template <ElementOrder O>
H5::DataSet SimulationAcFemFreqD3<O>::Impl::_begin_hdf5_file(
) {
    // create the file
    _hdf5_file = H5::H5File(_hdf5_file_path, H5F_ACC_TRUNC);

    // write signature
    write_string_attr(_hdf5_file, "Conventions", HDF5_SIGNATURE);

    // write simulation type
    H5::Group sim_type_grp = _hdf5_file.createGroup("/simulation_type");
    write_string_attr(sim_type_grp,
        "phenomenon", HDF5_PHENOMENON<Phenomenon::ACOUSTIC>
    );
    write_string_attr(sim_type_grp, 
        "numerical_method", HDF5_NUMERICAL_METHOD<NumericalMethod::FEM>
    );
    write_string_attr(sim_type_grp,
        "domain", HDF5_DOMAIN<Domain::FREQUENCY>
    );
    write_string_attr(sim_type_grp,
        "dimension", HDF5_DIMENSION<Dimension::D3>
    );
    write_string_attr(sim_type_grp,
        "element_order", HDF5_ORDER<O>
    );

    // allocate pressure solution matrix
    H5::Group solution_grp = _hdf5_file.createGroup("/results");
    H5::DataSet pressure_data_set = create_cmplx128_dataset(
        solution_grp,
        "pressure",
        _freq_count,
        _ni_count
    );
    write_string_attr(pressure_data_set, "units", "Pa");
    H5DSset_label(pressure_data_set.getId(), 0, "frequency_index");
    H5DSset_label(pressure_data_set.getId(), 1, "node_index");

    return pressure_data_set;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_write_simulation_inputs_to_hdf5_file(
) {
    // write simulated frequncies
    {
        H5::DataSet data_set = write_dataset_float64_1d(
            _hdf5_file,
            "simulated_frequencies",
            static_cast_contiguous_data<double>(
                _freq_steps.data(), _freq_count
            ).get(),
            _freq_count
        );
        write_string_attr(data_set, "units", "Hz");
    }

    // write mesh
    H5::Group mesh_grp = _hdf5_file.createGroup("/mesh");
    { // nodes
        H5::DataSet data_set = write_dataset_float64_2d(
            mesh_grp,
            "nodes",
            static_cast_contiguous_data<double>(
                _ni_to_coords.data()->data(),
                _ni_count * 3UL
            ).get(),
            _ni_count,
            3UL
        );
        write_string_attr(data_set, "units", "m");
    }
    { // surface_elements
        H5::DataSet data_set = write_dataset_uint64_2d(
            mesh_grp,
            "surface_elements",
            to_one_based_index(
                _sei_to_ni.data()->data(),
                _sei_count * NODES_IN_SFC_ELEM<O>
            ).get(),
            _sei_count,
            NODES_IN_SFC_ELEM<O>
        );
        write_string_attr(data_set, "units", "dimensionless");
        write_int_attr(data_set, "node_index_base", 1UL);
    }
    { // volume_elements
        H5::DataSet data_set = write_dataset_uint64_2d(
            mesh_grp,
            "volume_elements",
            to_one_based_index(
                _vei_to_ni.data()->data(),
                _vei_count * NODES_IN_VOL_ELEM<O>
            ).get(),
            _vei_count,
            NODES_IN_VOL_ELEM<O>
        );
        write_string_attr(data_set, "units", "dimensionless");
        write_int_attr(data_set, "node_index_base", 1UL);
    }
}

template<ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_write_solution_for_one_freq(
    H5::DataSet& pressure_data_set,
    const size_t fi
) {
    H5::CompType cmplx_type = get_hdf5_cmplx_type();
    const hsize_t row_size[2UL] = { 1UL, _ni_count };

    const hsize_t offset[2UL] = { fi, 0UL };
    H5::DataSpace matrix_data_space = pressure_data_set.getSpace();
    matrix_data_space.selectHyperslab(H5S_SELECT_SET, row_size, offset);
    
    const H5::DataSpace x_data_space(2UL, row_size);

    pressure_data_set.write(
        _x.data(), cmplx_type, x_data_space, matrix_data_space
    );
}

} // namespace numav

NUMAV_INSTANTIATE_SIM_AC_FEM_FREQ_D3
