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
    H5::Group inputs_grp = _hdf5_file.createGroup("/inputs");

    // simulated frequncies
    {
        H5::DataSet data_set = write_dataset_float64_1d(
            inputs_grp,
            "simulated_frequencies",
            static_cast_contiguous_data<double>(
                _freq_steps.data(), _freq_count
            ).get(),
            _freq_count
        );
        write_string_attr(data_set, "units", "Hz");
        H5DSset_label(data_set.getId(), 0, "frequency_index");
    }

    // mesh
    H5::Group mesh_grp = _hdf5_file.createGroup("/inputs/mesh");
    {   // nodes
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
        H5DSset_label(data_set.getId(), 0, "node_index");
        H5DSset_label(data_set.getId(), 1, "coordinates");
    }
    {   // surface_elements
        H5::DataSet data_set = write_dataset_uint64_2d(
            mesh_grp,
            "surface_elements",
            to_one_based_index(
                _sei_to_ni.data()->data(),
                _sei_count * ENI_COUNT_SFC<O>
            ).get(),
            _sei_count,
            ENI_COUNT_SFC<O>
        );
        write_string_attr(data_set, "units", "dimensionless");
        write_int_attr(data_set, "node_index_base", 1UL);
        H5DSset_label(data_set.getId(), 0, "surface_element_index");
        H5DSset_label(data_set.getId(), 1, "elementary_node_index");
    }
    {   // volume_elements
        H5::DataSet data_set = write_dataset_uint64_2d(
            mesh_grp,
            "volume_elements",
            to_one_based_index(
                _vei_to_ni.data()->data(),
                _vei_count * ENI_COUNT_VOL<O>
            ).get(),
            _vei_count,
            ENI_COUNT_VOL<O>
        );
        write_string_attr(data_set, "units", "dimensionless");
        write_int_attr(data_set, "node_index_base", 1UL);
        H5DSset_label(data_set.getId(), 0, "volume_element_index");
        H5DSset_label(data_set.getId(), 1, "elementary_node_index");
    }
    {   // surface_elements_material
        H5::DataSet data_set = write_dataset_uint64_1d(
            mesh_grp,
            "surface_elements_materials",
            _sei_to_espg.data(),
            _sei_count
        );
        write_string_attr(data_set, "units", "dimensionless");
        H5DSset_label(data_set.getId(), 0, "surface_element_index");
    }
    {   // volume_elements_material
        H5::DataSet data_set = write_dataset_uint64_1d(
            mesh_grp,
            "volume_elements_materials",
            _vei_to_evpg.data(),
            _vei_count
        );
        write_string_attr(data_set, "units", "dimensionless");
        H5DSset_label(data_set.getId(), 0, "volume_element_index");
    }

    // volume_materials
    H5::Group vol_mat_grp = _hdf5_file.createGroup("/inputs/volume_materials");
    {   // physical_groups
        fz::SafePtr<uint64_t> _ivpg_to_evpg(_ivpg_count);
        for (size_t ivpg = 0UL; ivpg != _ivpg_count; ++ivpg) {
            const size_t evpg = _evpg_ivpg_bimap.right.at(ivpg);
            _ivpg_to_evpg[ivpg] = evpg;
        }
        H5::DataSet data_set = write_dataset_uint64_1d(
            vol_mat_grp,
            "physical_groups",
            _ivpg_to_evpg.data(),
            _ivpg_count
        );
        _ivpg_to_evpg.free();
        write_string_attr(data_set, "units", "dimensionless");
        H5DSset_label(data_set.getId(), 0, "volume_material_index");
    }
    {   // density
        fz::SafePtr<Cmplx> density(_ivpg_count * _freq_count);
        for (size_t ivpg = 0UL; ivpg != _ivpg_count; ++ivpg) {
            for (size_t fi = 0UL; fi != _freq_count; ++fi) {
                const Float freq = _freq_steps[fi];
                density[ivpg*_freq_count + fi] = 
                    (_ivpg_to_volprop[ivpg].density)(freq);
            }
        }
        H5::DataSet data_set = write_dataset_cmplx128_2d(
            vol_mat_grp,
            "density",
            density.data(),
            _ivpg_count,
            _freq_count
        );
        density.free();
        write_string_attr(data_set, "units", "kg/m^3");
        H5DSset_label(data_set.getId(), 0, "volume_material_index");
        H5DSset_label(data_set.getId(), 1, "frequency_index");
    }
    {   // sound_speed
        fz::SafePtr<Cmplx> soundspeed(_ivpg_count * _freq_count);
        for (size_t ivpg = 0UL; ivpg != _ivpg_count; ++ivpg) {
            for (size_t fi = 0UL; fi != _freq_count; ++fi) {
                const Float freq = _freq_steps[fi];
                soundspeed[ivpg*_freq_count + fi] = 
                    (_ivpg_to_volprop[ivpg].soundspeed)(freq);
            }
        }
        H5::DataSet data_set = write_dataset_cmplx128_2d(
            vol_mat_grp,
            "sound_speed",
            soundspeed.data(),
            _ivpg_count,
            _freq_count
        );
        soundspeed.free();
        write_string_attr(data_set, "units", "kg/m^3");
        H5DSset_label(data_set.getId(), 0, "volume_material_index");
        H5DSset_label(data_set.getId(), 1, "frequency_index");
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
