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
        _fi_count,
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
                _fi_to_freq.data(), _fi_count
            ).get(),
            _fi_count
        );
        write_string_attr(data_set, "units", "Hz");
        H5DSset_label(data_set.getId(), 0, "frequency_index");
    }

    // mesh
    H5::Group mesh_grp = inputs_grp.createGroup("mesh");
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
                _sei_count * ENIS_COUNT<O>
            ).get(),
            _sei_count,
            ENIS_COUNT<O>
        );
        write_string_attr(data_set, "units", "dimensionless");
        H5DSset_label(data_set.getId(), 0, "surface_element_index");
        H5DSset_label(data_set.getId(), 1, "elemental_node_index");
        write_int_attr(data_set, "node_index_base", 1UL);
    }
    {   // volume_elements
        H5::DataSet data_set = write_dataset_uint64_2d(
            mesh_grp,
            "volume_elements",
            to_one_based_index(
                _vei_to_ni.data()->data(),
                _vei_count * ENIV_COUNT<O>
            ).get(),
            _vei_count,
            ENIV_COUNT<O>
        );
        write_string_attr(data_set, "units", "dimensionless");
        H5DSset_label(data_set.getId(), 0, "volume_element_index");
        H5DSset_label(data_set.getId(), 1, "elemental_node_index");
        write_int_attr(data_set, "node_index_base", 1UL);
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
    H5::Group vol_mat_grp = inputs_grp.createGroup("volume_materials");
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
        fz::SafePtr<Cmplx> density(_ivpg_count * _fi_count);
        for (size_t ivpg = 0UL; ivpg != _ivpg_count; ++ivpg) {
            for (size_t fi = 0UL; fi != _fi_count; ++fi) {
                const Float freq = _fi_to_freq[fi];
                density[ivpg*_fi_count + fi] = 
                    (_ivpg_to_volprop[ivpg].density)(freq);
            }
        }
        H5::DataSet data_set = write_dataset_cmplx128_2d(
            vol_mat_grp,
            "density",
            density.data(),
            _ivpg_count,
            _fi_count
        );
        density.free();
        write_string_attr(data_set, "units", "kg/m3");
        H5DSset_label(data_set.getId(), 0, "volume_material_index");
        H5DSset_label(data_set.getId(), 1, "frequency_index");
    }
    {   // sound_speed
        fz::SafePtr<Cmplx> soundspeed(_ivpg_count * _fi_count);
        for (size_t ivpg = 0UL; ivpg != _ivpg_count; ++ivpg) {
            for (size_t fi = 0UL; fi != _fi_count; ++fi) {
                const Float freq = _fi_to_freq[fi];
                soundspeed[ivpg*_fi_count + fi] = 
                    (_ivpg_to_volprop[ivpg].soundspeed)(freq);
            }
        }
        H5::DataSet data_set = write_dataset_cmplx128_2d(
            vol_mat_grp,
            "sound_speed",
            soundspeed.data(),
            _ivpg_count,
            _fi_count
        );
        soundspeed.free();
        write_string_attr(data_set, "units", "m/s");
        H5DSset_label(data_set.getId(), 0, "volume_material_index");
        H5DSset_label(data_set.getId(), 1, "frequency_index");
    }

    // surface_materials
    H5::Group sfc_mat_grp = inputs_grp.createGroup("surface_materials");
    {   // physical_groups
        fz::SafePtr<uint64_t> _ispgi_to_espg(_ispgi_count);
        for (size_t ispgi = 0UL; ispgi != _ispgi_count; ++ispgi) {
            const size_t espg = _espg_ispgi_bimap.right.at(ispgi);
            _ispgi_to_espg[ispgi] = espg;
        }
        H5::DataSet data_set = write_dataset_uint64_1d(
            sfc_mat_grp,
            "physical_groups",
            _ispgi_to_espg.data(),
            _ispgi_count
        );
        _ispgi_to_espg.free();
        write_string_attr(data_set, "units", "dimensionless");
        H5DSset_label(data_set.getId(), 0, "surface_material_index");
    }
    {   // impedance
        fz::SafePtr<Cmplx> impedance(_ispgi_count * _fi_count);
        for (size_t ispgi = 0UL; ispgi != _ispgi_count; ++ispgi) {
            for (size_t fi = 0UL; fi != _fi_count; ++fi) {
                const Float freq = _fi_to_freq[fi];
                impedance[ispgi*_fi_count + fi] = 
                    (_ispgi_to_impedance[ispgi])(freq);
            }
        }
        H5::DataSet data_set = write_dataset_cmplx128_2d(
            sfc_mat_grp,
            "impedance",
            impedance.data(),
            _ispgi_count,
            _fi_count
        );
        impedance.free();
        write_string_attr(data_set, "units", "Pa.s/m");
        H5DSset_label(data_set.getId(), 0, "surface_material_index");
        H5DSset_label(data_set.getId(), 1, "frequency_index");
    }

    // sound_sources
    H5::Group src_grp = inputs_grp.createGroup("sound_sources");
    H5::Group pnt_src_grp = src_grp.createGroup("point");
    H5::Group pnt_volvel_grp = pnt_src_grp.createGroup("volume_velocity");
    {   // point node_index
        H5::DataSet data_set = write_dataset_uint64_1d(
            pnt_volvel_grp,
            "node_index",
            to_one_based_index(_vpi_to_ni.data(), _vpi_count).get(),
            _vpi_count
        );
        write_string_attr(data_set, "units", "m");
        H5DSset_label(data_set.getId(), 0, "volume_velocity_point_index");
        write_int_attr(data_set, "node_index_base", 1UL);
    }
    {   // point volume_velocity 
        fz::SafePtr<Cmplx> volvel(_vpi_count * _fi_count);
        for (size_t vpi = 0UL; vpi != _vpi_count; ++vpi) {
            for (size_t fi = 0UL; fi != _fi_count; ++fi) {
                const Float freq = _fi_to_freq[fi];
                volvel[vpi*_fi_count + fi] = (_vpi_to_volvel[vpi])(freq);
            }
        }
        H5::DataSet data_set = write_dataset_cmplx128_2d(
            pnt_volvel_grp,
            "volume_velocity",
            volvel.data(),
            _vpi_count,
            _fi_count
        );
        volvel.free();
        write_string_attr(data_set, "units", "m3/s");
        H5DSset_label(data_set.getId(), 0, "volume_velocity_point_index");
        H5DSset_label(data_set.getId(), 1, "frequency_index");
    }
    H5::Group pnt_pressure_grp = pnt_src_grp.createGroup("pressure");
    {   // point node_index
        H5::DataSet data_set = write_dataset_uint64_1d(
            pnt_pressure_grp,
            "node_index",
            to_one_based_index(_ppi_to_ni.data(), _ppi_count).get(),
            _ppi_count
        );
        write_string_attr(data_set, "units", "m");
        H5DSset_label(data_set.getId(), 0, "pressure_point_index");
        write_int_attr(data_set, "node_index_base", 1UL);
    }
    {   // point pressure
        fz::SafePtr<Cmplx> pressure(_ppi_count * _fi_count);
        for (size_t ppi = 0UL; ppi != _ppi_count; ++ppi) {
            for (size_t fi = 0UL; fi != _fi_count; ++fi) {
                const Float freq = _fi_to_freq[fi];
                pressure[ppi*_fi_count + fi] = (_ppi_to_pressure[ppi])(freq);
            }
        }
        H5::DataSet data_set = write_dataset_cmplx128_2d(
            pnt_pressure_grp,
            "pressure",
            pressure.data(),
            _ppi_count,
            _fi_count
        );
        pressure.free();
        write_string_attr(data_set, "units", "Pa");
        H5DSset_label(data_set.getId(), 0, "pressure_point_index");
        H5DSset_label(data_set.getId(), 1, "frequency_index");
    }

    H5::Group sfc_src_grp = src_grp.createGroup("surface");
    H5::Group sfc_velocity_grp = sfc_src_grp.createGroup("particle_velocity");
    {   // surface physical_groups
        fz::SafePtr<uint64_t> ispgv_to_espg(_ispgv_count);
        for (size_t ispgv = 0UL; ispgv != _ispgv_count; ++ispgv) {
            ispgv_to_espg[ispgv] = _espg_ispgv_bimap.right.at(ispgv);
        }
        H5::DataSet data_set = write_dataset_uint64_1d(
            sfc_velocity_grp,
            "physical_groups",
            ispgv_to_espg.data(),
            _ispgv_count
        );
        ispgv_to_espg.free();
        write_string_attr(data_set, "units", "dimensionless");
        H5DSset_label(data_set.getId(), 0, "particle_velocity_surface_index");
    }
    {   // surface particle_velocity 
        fz::SafePtr<Cmplx> vel(_ispgv_count * _fi_count);
        for (size_t ispgv = 0UL; ispgv != _ispgv_count; ++ispgv) {
            for (size_t fi = 0UL; fi != _fi_count; ++fi) {
                const Float freq = _fi_to_freq[fi];
                vel[ispgv*_fi_count + fi] = (_ispgv_to_velocity[ispgv])(freq);
            }
        }
        H5::DataSet data_set = write_dataset_cmplx128_2d(
            sfc_velocity_grp,
            "particle_velocity",
            vel.data(),
            _ispgv_count,
            _fi_count
        );
        vel.free();
        write_string_attr(data_set, "units", "m/s");
        H5DSset_label(data_set.getId(), 0, "particle_velocity_surface_index");
        H5DSset_label(data_set.getId(), 1, "frequency_index");
    }
    H5::Group sfc_pressure_grp = sfc_src_grp.createGroup("pressure");
    {   // surface physical_groups
        fz::SafePtr<uint64_t> ispgp_to_espg(_ispgp_count);
        for (size_t ispgp = 0UL; ispgp != _ispgp_count; ++ispgp) {
            ispgp_to_espg[ispgp] = _espg_ispgp_bimap.right.at(ispgp);
        }
        H5::DataSet data_set = write_dataset_uint64_1d(
            sfc_pressure_grp,
            "physical_groups",
            ispgp_to_espg.data(),
            _ispgp_count
        );
        ispgp_to_espg.free();
        write_string_attr(data_set, "units", "dimensionless");
        H5DSset_label(data_set.getId(), 0, "pressure_surface_index");
    }
    {   // surface pressure 
        fz::SafePtr<Cmplx> pres(_ispgp_count * _fi_count);
        for (size_t ispgp = 0UL; ispgp != _ispgp_count; ++ispgp) {
            for (size_t fi = 0UL; fi != _fi_count; ++fi) {
                const Float freq = _fi_to_freq[fi];
                pres[ispgp*_fi_count + fi] = (_ispgp_to_pressure[ispgp])(freq);
            }
        }
        H5::DataSet data_set = write_dataset_cmplx128_2d(
            sfc_pressure_grp,
            "pressure",
            pres.data(),
            _ispgp_count,
            _fi_count
        );
        pres.free();
        write_string_attr(data_set, "units", "m/s");
        H5DSset_label(data_set.getId(), 0, "pressure_surface_index");
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
