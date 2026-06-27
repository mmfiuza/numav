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
H5::DataSet SimulationAcFemFreqD3Tet<O>::Impl::_begin_hdf5_file(
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
    H5::DataSet pressure_dataset = create_cmplx128_dataset(
        solution_grp,
        "pressure",
        _fi_count,
        _ni_count
    );
    write_string_attr(pressure_dataset, "units", "Pa");
    H5DSset_label(pressure_dataset.getId(), 0UL, "frequency_index");
    H5DSset_label(pressure_dataset.getId(), 1UL, "node_index");

    return pressure_dataset;
}

template <ElementOrder O>
void SimulationAcFemFreqD3Tet<O>::Impl::_write_simulation_inputs_to_hdf5_file(
) {
    H5::Group inputs_grp = _hdf5_file.createGroup("/inputs");

    // simulated frequencies
    H5::DataSet dataset_freq = write_dataset_float64_1d(
        inputs_grp,
        "simulated_frequencies",
        static_cast_contiguous_data<double>(
            _fi_to_freq.data(), _fi_to_freq.size()
        ).get(),
        _fi_count
    );
    write_string_attr(dataset_freq, "units", "Hz");
    H5DSset_label(dataset_freq.getId(), 0UL, "frequency_index");
    H5DSset_scale(dataset_freq.getId(), "simulated_frequencies");

    // mesh
    H5::Group mesh_grp = inputs_grp.createGroup("mesh");
    {   // nodes
        H5::DataSet dataset = write_dataset_float64_2d(
            mesh_grp,
            "nodes",
            static_cast_contiguous_data<double>(
                _ni_to_coords.data()->data(),
                _ni_count * 3UL
            ).get(),
            _ni_count,
            3UL
        );
        write_string_attr(dataset, "units", "m");
        H5DSset_label(dataset.getId(), 0UL, "node_index");
        H5DSset_label(dataset.getId(), 1UL, "coordinates");
    }
    {   // surface_elements
        H5::DataSet dataset = write_dataset_uint64_2d(
            mesh_grp,
            "surface_elements",
            to_one_based_index(
                _sei_to_ni.data()->data(),
                _sei_count * ENIS_COUNT<O>
            ).get(),
            _sei_count,
            ENIS_COUNT<O>
        );
        write_string_attr(dataset, "units", "1");
        H5DSset_label(dataset.getId(), 0UL, "surface_element_index");
        H5DSset_label(dataset.getId(), 1UL, "elemental_node_index");
        write_uint64_attr(dataset, "node_index_base", 1UL);
    }
    {   // volume_elements
        H5::DataSet dataset = write_dataset_uint64_2d(
            mesh_grp,
            "volume_elements",
            to_one_based_index(
                _vei_to_ni.data()->data(),
                _vei_count * ENIV_COUNT<O>
            ).get(),
            _vei_count,
            ENIV_COUNT<O>
        );
        write_string_attr(dataset, "units", "1");
        H5DSset_label(dataset.getId(), 0UL, "volume_element_index");
        H5DSset_label(dataset.getId(), 1UL, "elemental_node_index");
        write_uint64_attr(dataset, "node_index_base", 1UL);
    }
    {   // surface_elements_material
        H5::DataSet dataset = write_dataset_uint64_1d(
            mesh_grp,
            "surface_elements_materials",
            _sei_to_espg.data(),
            _sei_count
        );
        write_string_attr(dataset, "units", "1");
        H5DSset_label(dataset.getId(), 0UL, "surface_element_index");
    }
    {   // volume_elements_material
        H5::DataSet dataset = write_dataset_uint64_1d(
            mesh_grp,
            "volume_elements_materials",
            _vei_to_evpg.data(),
            _vei_count
        );
        write_string_attr(dataset, "units", "1");
        H5DSset_label(dataset.getId(), 0UL, "volume_element_index");
    }

    // volume_materials
    H5::Group vol_mat_grp = inputs_grp.createGroup("volume_materials");
    {   // physical_groups
        fz::SafePtr<uint64_t> _ivpg_to_evpg(_ivpg_count);
        for (uint64_t ivpg = 0UL; ivpg != _ivpg_count; ++ivpg) {
            const uint64_t evpg = _evpg_ivpg_bimap.right.at(ivpg);
            _ivpg_to_evpg[ivpg] = evpg;
        }
        H5::DataSet dataset = write_dataset_uint64_1d(
            vol_mat_grp,
            "physical_groups",
            _ivpg_to_evpg.data(),
            _ivpg_count
        );
        _ivpg_to_evpg.free();
        write_string_attr(dataset, "units", "1");
        H5DSset_label(dataset.getId(), 0UL, "volume_material_index");
    }
    {   // density
        fz::SafePtr<Cmplx> density(_ivpg_count * _fi_count);
        for (uint64_t ivpg = 0UL; ivpg != _ivpg_count; ++ivpg) {
            for (uint64_t fi = 0UL; fi != _fi_count; ++fi) {
                const Float freq = _fi_to_freq[fi];
                density[ivpg*_fi_count + fi] = 
                    (_ivpg_to_volprop[ivpg].density)(freq);
            }
        }
        H5::DataSet dataset = write_dataset_cmplx128_2d(
            vol_mat_grp,
            "density",
            static_cast_contiguous_data<std::complex<double>>(
                density.data(), density.size()
            ).get(),
            _ivpg_count,
            _fi_count
        );
        density.free();
        write_string_attr(dataset, "units", "kg m-3");
        H5DSset_label(dataset.getId(), 0UL, "volume_material_index");
        H5DSset_label(dataset.getId(), 1UL, "frequency_index");
        H5DSattach_scale(dataset.getId(), dataset_freq.getId(), 1UL);
    }
    {   // sound_speed
        fz::SafePtr<Cmplx> soundspeed(_ivpg_count * _fi_count);
        for (uint64_t ivpg = 0UL; ivpg != _ivpg_count; ++ivpg) {
            for (uint64_t fi = 0UL; fi != _fi_count; ++fi) {
                const Float freq = _fi_to_freq[fi];
                soundspeed[ivpg*_fi_count + fi] = 
                    (_ivpg_to_volprop[ivpg].soundspeed)(freq);
            }
        }
        H5::DataSet dataset = write_dataset_cmplx128_2d(
            vol_mat_grp,
            "sound_speed",
            static_cast_contiguous_data<std::complex<double>>(
                soundspeed.data(), soundspeed.size()
            ).get(),
            _ivpg_count,
            _fi_count
        );
        soundspeed.free();
        write_string_attr(dataset, "units", "m/s");
        H5DSset_label(dataset.getId(), 0UL, "volume_material_index");
        H5DSset_label(dataset.getId(), 1UL, "frequency_index");
        H5DSattach_scale(dataset.getId(), dataset_freq.getId(), 1UL);
    }

    // surface_materials
    H5::Group sfc_mat_grp = inputs_grp.createGroup("surface_materials");
    {   // physical_groups
        fz::SafePtr<uint64_t> _ispgi_to_espg(_ispgi_count);
        for (uint64_t ispgi = 0UL; ispgi != _ispgi_count; ++ispgi) {
            const uint64_t espg = _espg_ispgi_bimap.right.at(ispgi);
            _ispgi_to_espg[ispgi] = espg;
        }
        H5::DataSet dataset = write_dataset_uint64_1d(
            sfc_mat_grp,
            "physical_groups",
            _ispgi_to_espg.data(),
            _ispgi_count
        );
        _ispgi_to_espg.free();
        write_string_attr(dataset, "units", "1");
        H5DSset_label(dataset.getId(), 0UL, "surface_material_index");
    }
    {   // impedance
        fz::SafePtr<Cmplx> impedance(_ispgi_count * _fi_count);
        for (uint64_t ispgi = 0UL; ispgi != _ispgi_count; ++ispgi) {
            for (uint64_t fi = 0UL; fi != _fi_count; ++fi) {
                const Float freq = _fi_to_freq[fi];
                impedance[ispgi*_fi_count + fi] = 
                    (_ispgi_to_impedance[ispgi])(freq);
            }
        }
        H5::DataSet dataset = write_dataset_cmplx128_2d(
            sfc_mat_grp,
            "impedance",
            static_cast_contiguous_data<std::complex<double>>(
                impedance.data(), impedance.size()
            ).get(),
            _ispgi_count,
            _fi_count
        );
        impedance.free();
        write_string_attr(dataset, "units", "Pa s m-1");
        H5DSset_label(dataset.getId(), 0UL, "surface_material_index");
        H5DSset_label(dataset.getId(), 1UL, "frequency_index");
        H5DSattach_scale(dataset.getId(), dataset_freq.getId(), 1UL);
    }

    // sound_sources
    H5::Group src_grp = inputs_grp.createGroup("sound_sources");
    H5::Group pnt_src_grp = src_grp.createGroup("point");
    H5::Group pnt_volvel_grp = pnt_src_grp.createGroup("volume_velocity");
    {   // point node_index
        H5::DataSet dataset = write_dataset_uint64_1d(
            pnt_volvel_grp,
            "node_index",
            to_one_based_index(_vpi_to_ni.data(), _vpi_count).get(),
            _vpi_count
        );
        write_string_attr(dataset, "units", "1");
        H5DSset_label(dataset.getId(), 0UL, "volume_velocity_point_index");
        write_uint64_attr(dataset, "node_index_base", 1UL);
    }
    {   // point volume_velocity 
        fz::SafePtr<Cmplx> volvel(_vpi_count * _fi_count);
        for (uint64_t vpi = 0UL; vpi != _vpi_count; ++vpi) {
            for (uint64_t fi = 0UL; fi != _fi_count; ++fi) {
                const Float freq = _fi_to_freq[fi];
                volvel[vpi*_fi_count + fi] = (_vpi_to_volvel[vpi])(freq);
            }
        }
        H5::DataSet dataset = write_dataset_cmplx128_2d(
            pnt_volvel_grp,
            "volume_velocity",
            static_cast_contiguous_data<std::complex<double>>(
                volvel.data(), volvel.size()
            ).get(),
            _vpi_count,
            _fi_count
        );
        volvel.free();
        write_string_attr(dataset, "units", "m3 s-1");
        H5DSset_label(dataset.getId(), 0UL, "volume_velocity_point_index");
        H5DSset_label(dataset.getId(), 1UL, "frequency_index");
        H5DSattach_scale(dataset.getId(), dataset_freq.getId(), 1UL);
    }
    H5::Group pnt_pressure_grp = pnt_src_grp.createGroup("pressure");
    {   // point node_index
        H5::DataSet dataset = write_dataset_uint64_1d(
            pnt_pressure_grp,
            "node_index",
            to_one_based_index(_ppi_to_ni.data(), _ppi_count).get(),
            _ppi_count
        );
        write_string_attr(dataset, "units", "1");
        H5DSset_label(dataset.getId(), 0UL, "pressure_point_index");
        write_uint64_attr(dataset, "node_index_base", 1UL);
    }
    {   // point pressure
        fz::SafePtr<Cmplx> pressure(_ppi_count * _fi_count);
        for (uint64_t ppi = 0UL; ppi != _ppi_count; ++ppi) {
            for (uint64_t fi = 0UL; fi != _fi_count; ++fi) {
                const Float freq = _fi_to_freq[fi];
                pressure[ppi*_fi_count + fi] = (_ppi_to_pressure[ppi])(freq);
            }
        }
        H5::DataSet dataset = write_dataset_cmplx128_2d(
            pnt_pressure_grp,
            "pressure",
            static_cast_contiguous_data<std::complex<double>>(
                pressure.data(), pressure.size()
            ).get(),
            _ppi_count,
            _fi_count
        );
        pressure.free();
        write_string_attr(dataset, "units", "Pa");
        H5DSset_label(dataset.getId(), 0UL, "pressure_point_index");
        H5DSset_label(dataset.getId(), 1UL, "frequency_index");
        H5DSattach_scale(dataset.getId(), dataset_freq.getId(), 1UL);
    }

    H5::Group sfc_src_grp = src_grp.createGroup("surface");
    H5::Group sfc_velocity_grp = sfc_src_grp.createGroup("particle_velocity");
    {   // surface physical_groups
        fz::SafePtr<uint64_t> ispgv_to_espg(_ispgv_count);
        for (uint64_t ispgv = 0UL; ispgv != _ispgv_count; ++ispgv) {
            ispgv_to_espg[ispgv] = _espg_ispgv_bimap.right.at(ispgv);
        }
        H5::DataSet dataset = write_dataset_uint64_1d(
            sfc_velocity_grp,
            "physical_groups",
            ispgv_to_espg.data(),
            _ispgv_count
        );
        ispgv_to_espg.free();
        write_string_attr(dataset, "units", "1");
        H5DSset_label(dataset.getId(), 0UL, "particle_velocity_surface_index");
    }
    {   // surface particle_velocity 
        fz::SafePtr<Cmplx> vel(_ispgv_count * _fi_count);
        for (uint64_t ispgv = 0UL; ispgv != _ispgv_count; ++ispgv) {
            for (uint64_t fi = 0UL; fi != _fi_count; ++fi) {
                const Float freq = _fi_to_freq[fi];
                vel[ispgv*_fi_count + fi] = (_ispgv_to_velocity[ispgv])(freq);
            }
        }
        H5::DataSet dataset = write_dataset_cmplx128_2d(
            sfc_velocity_grp,
            "particle_velocity",
            static_cast_contiguous_data<std::complex<double>>(
                vel.data(), vel.size()
            ).get(),
            _ispgv_count,
            _fi_count
        );
        vel.free();
        write_string_attr(dataset, "units", "m/s");
        H5DSset_label(dataset.getId(), 0UL, "particle_velocity_surface_index");
        H5DSset_label(dataset.getId(), 1UL, "frequency_index");
        H5DSattach_scale(dataset.getId(), dataset_freq.getId(), 1UL);
    }
    H5::Group sfc_pressure_grp = sfc_src_grp.createGroup("pressure");
    {   // surface physical_groups
        fz::SafePtr<uint64_t> ispgp_to_espg(_ispgp_count);
        for (uint64_t ispgp = 0UL; ispgp != _ispgp_count; ++ispgp) {
            ispgp_to_espg[ispgp] = _espg_ispgp_bimap.right.at(ispgp);
        }
        H5::DataSet dataset = write_dataset_uint64_1d(
            sfc_pressure_grp,
            "physical_groups",
            ispgp_to_espg.data(),
            _ispgp_count
        );
        ispgp_to_espg.free();
        write_string_attr(dataset, "units", "1");
        H5DSset_label(dataset.getId(), 0UL, "pressure_surface_index");
    }
    {   // surface pressure 
        fz::SafePtr<Cmplx> pres(_ispgp_count * _fi_count);
        for (uint64_t ispgp = 0UL; ispgp != _ispgp_count; ++ispgp) {
            for (uint64_t fi = 0UL; fi != _fi_count; ++fi) {
                const Float freq = _fi_to_freq[fi];
                pres[ispgp*_fi_count + fi] = (_ispgp_to_pressure[ispgp])(freq);
            }
        }
        H5::DataSet dataset = write_dataset_cmplx128_2d(
            sfc_pressure_grp,
            "pressure",
            static_cast_contiguous_data<std::complex<double>>(
                pres.data(), pres.size()
            ).get(),
            _ispgp_count,
            _fi_count
        );
        pres.free();
        write_string_attr(dataset, "units", "Pa");
        H5DSset_label(dataset.getId(), 0UL, "pressure_surface_index");
        H5DSset_label(dataset.getId(), 1UL, "frequency_index");
        H5DSattach_scale(dataset.getId(), dataset_freq.getId(), 1UL);
    }
}

template<ElementOrder O>
void SimulationAcFemFreqD3Tet<O>::Impl::_write_solution_for_one_freq(
    H5::DataSet& pressure_dataset,
    const uint64_t fi
) {
    H5::CompType cmplx_type = get_hdf5_cmplx_type();
    const hsize_t row_size[2UL] = { 1UL, _ni_count };

    const hsize_t offset[2UL] = { fi, 0UL };
    H5::DataSpace matrix_data_space = pressure_dataset.getSpace();
    matrix_data_space.selectHyperslab(H5S_SELECT_SET, row_size, offset);
    
    const H5::DataSpace x_data_space(2UL, row_size);

    pressure_dataset.write(
        static_cast_contiguous_data<std::complex<double>>(
            _x.data(), _x.size()
        ).get(),
        cmplx_type,
        x_data_space,
        matrix_data_space
    );
}

} // namespace numav

NUMAV_INSTANTIATE_SIM_AC_FEM_FREQ_D3
