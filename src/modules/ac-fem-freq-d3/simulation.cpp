// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include <tuple>
#include <cmath>
#include <charconv>
#include <fstream>
#undef NDEBUG
#include <cassert>

#include "Eigen/Eigen"

#include "numav/numav.hpp"
#include "numav/typedefs.hpp"
#include "common/log.hpp"
#include "common/hash_functions.hpp"
#include "common/maths.hpp"

namespace numav {

static constexpr size_t DIM = DIM_COUNT<Dimension::D3>;
static constexpr double VOLUME_REF_TET = 1/6;

using ResultAcFemFreqD3 = typename numav::Result<
    Phenomenon::ACOUSTIC,
    NumericalMethod::FEM,
    Domain::FREQUENCY,
    Dimension::D3
>;

ResultAcFemFreqD3::Result() = default;
ResultAcFemFreqD3::~Result() = default;

template<ElementOrder O>
using SimulationAcFemFreqD3 = typename numav::Simulation<
    Phenomenon::ACOUSTIC,
    NumericalMethod::FEM,
    Domain::FREQUENCY,
    Dimension::D3,
    O
>;

template<ElementOrder O>
SimulationAcFemFreqD3<O>::Simulation() {
    log::set_level();
    log::set_pattern();
    _is_mesh_defined = false;
    _is_freq_range_defined = false;
    _is_any_source_defined = false;
}

template<ElementOrder O>
SimulationAcFemFreqD3<O>::~Simulation() {
    _freq_vec.free();
    _node_coords.free();
    _sfc_elem_node_idx.free();
    _vol_elem_node_idx.free();
    _epg_sfc_elem.free();
    _epg_vol_elem.free();
    // _ipg_sfc_elem.free();
    _ipg_vol_elem.free();
    _ipg_to_sfc_volvel.free();
    _ipg_to_sfc_pressure.free();
    _ipg_to_sfc_impedance.free();
    _ipg_to_volprop.free();
    _nonzero_row_idx.free();
    _nonzero_col_idx.free();
    _btb_detj_w_vals.free();
    // _nnt_detj_w_vals.free();
}

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::_node_count() const {
    return _node_coords.size();
}

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::_sfc_elem_count() const {
    return _sfc_elem_node_idx.size();
}

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::_vol_elem_count() const {
    return _vol_elem_node_idx.size();
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::_check_if_mesh_is_defined() {
    if (!_is_mesh_defined){
        log::error("Mesh not defined. Call load_mesh to do so.");
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::set_frequency_range(
    const double& freq_min,
    const double& freq_max
) {
    if (freq_min<0 || freq_max<0) {
        log::error("Frequency limits should be positive.");
    }
    if (freq_min >= freq_max) {
        log::error("Upper frequency should be greater than the lower.");
    }
    _freq_min = freq_min;
    _freq_max = freq_max;
    _is_freq_range_defined = true;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::load_mesh(
    const char* const path_to_mesh
) {
    if (_is_mesh_defined) {
        log::error("Mesh is already defined.");
    }
    std::string str = path_to_mesh;
    if (str.ends_with(".bdf") || str.ends_with(".nas")) {
        _load_bdf(path_to_mesh);
    }
    else {
        const size_t dot_position = str.find_last_of('.');
        const size_t format_len = str.size() - dot_position;
        const std::string format = std::string(
            str.substr(dot_position, format_len)
        );
        log::error("Unrecognized file format: \"{0}\".", format);
    }
    _generate_extra_nodes(); // call is based on the element order
    _is_mesh_defined = true;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_volume_material(
    const _epg_t& epg,
    const _FuncRealToCmplx& density,
    const _FuncRealToCmplx& soundspeed
) {
    _check_if_mesh_is_defined();
    if (!_existing_epg_vol.contains(epg)) {
        log::error("Tag {} not found in mesh file.", epg);
    }
    if (_epg_to_volprop.contains(epg)) {
        log::error("Tag {} already assigned.", epg);
    }
    _epg_to_volprop.insert({epg, {density, soundspeed}});
    const _ipg_t ipg = _epg_to_ipg_vol.size();
    _epg_to_ipg_vol.insert({epg, ipg});
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_sound_source(
    const TypeOfSource& type_of_source,
    const std::array<double,3>& point_coordinates,
    const PhysicalQuantity& physical_quantity_type,
    const _FuncRealToCmplx& physical_quantity_value
) {
    _check_if_mesh_is_defined();
    if (type_of_source != TypeOfSource::POINT) {
        log::error("Tried to assign coordinates to a surface sound source.");
    }
    const _idx_t closest_point_idx = _get_closest_point(point_coordinates);
    
    if (physical_quantity_type == PhysicalQuantity::VOLUME_VELOCITY) {
        _point_volvel.push_back(
            std::make_tuple(closest_point_idx, physical_quantity_value)
        );
    }
    else if (physical_quantity_type == PhysicalQuantity::PRESSURE) {
        _point_volvel.push_back(
            std::make_tuple(closest_point_idx, physical_quantity_value)
        );
    }
    else {
        log::error("Possible physical quantities are volume velocity"
            " or pressure.");
    }
    _is_any_source_defined = true;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_sound_source(
    const TypeOfSource& type_of_source,
    const _epg_t& epg,
    const PhysicalQuantity& physical_quantity_type,
    const _FuncRealToCmplx& physical_quantity_value
) {
    _check_if_mesh_is_defined();
    if (type_of_source != TypeOfSource::SURFACE) {
        log::error("Tried to assign a tag to a point.");
    }
    if (!_existing_epg_sfc.contains(epg)) {
        log::error("Tag {} not found in mesh file.", epg);
    }
    if (_epg_to_sfc_pressure.contains(epg) ||
        _epg_to_sfc_volvel.contains(epg) ||
        _epg_to_sfc_impedance.contains(epg)
    ) {
        log::error("Tag {} already assigned.", epg);
    }
    if (physical_quantity_type == PhysicalQuantity::PRESSURE) {
        const _ipg_t ipg = _epg_to_sfc_pressure.size();
        _epg_to_sfc_pressure.insert({epg, physical_quantity_value});
        _epg_to_ipg_sfc.insert({epg, ipg});
    }
    else if (physical_quantity_type == PhysicalQuantity::VOLUME_VELOCITY) {
        const _ipg_t ipg = _epg_to_sfc_volvel.size();
        _epg_to_sfc_volvel.insert({epg, physical_quantity_value});
        _epg_to_ipg_sfc.insert({epg, ipg});
    }
    else {
        log::error("Possible physical quantities are volume velocity"
            " or pressure.");
    }
    _is_any_source_defined = true;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::add_surface_specific_acoustic_impedance(
    const _epg_t& epg,
    const _FuncRealToCmplx& impedance
) {
    _check_if_mesh_is_defined();
    if (!_existing_epg_sfc.contains(epg)) {
        log::error("Tag {} not found in mesh file.", epg);
    }
    if (_epg_to_sfc_pressure.contains(epg) ||
        _epg_to_sfc_volvel.contains(epg) ||
        _epg_to_sfc_impedance.contains(epg)
    ) {
        log::error("Tag {} already assigned.", epg);
    }
    const _ipg_t ipg = _epg_to_sfc_impedance.size();
    _epg_to_sfc_impedance.insert({epg, impedance});
    _epg_to_ipg_sfc.insert({epg, ipg});
}

template <ElementOrder O>
_idx_t SimulationAcFemFreqD3<O>::_get_closest_point(
    const std::array<double,3>&
) {
    return 0; // TODO
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::_check_if_it_can_run() {
    _check_if_mesh_is_defined();
    if (!_is_any_source_defined){
        log::error("No sound source was defined."
            " Call add_sound_source to do so.");
    }
    if (!_is_freq_range_defined){
        log::error("Simulation frequency range was not defined."
            " Call set_frequency_range to do so.");
    }
    for (auto& epg : _existing_epg_vol) {
        if (!_epg_to_volprop.contains(epg)) {
            log::error("Volume tag {} was not assigned."
            " Call add_volume_material to do so.", epg);
        }
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::_define_freq_vector() {
    // todo: decide number here
    // todo: make it not linear
    _freq_vec = linspace(_freq_min, _freq_max, 8192);
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::_organize_physical_group_data()
{
    // generate structures accessed through IPG
    _ipg_to_sfc_volvel =
        fz::SafePtr<_FuncRealToCmplx>(_epg_to_sfc_volvel.size());
    for (const auto& [epg, volvel] : _epg_to_sfc_volvel) {
        const _ipg_t ipg = _epg_to_ipg_sfc.at(epg);
        _ipg_to_sfc_volvel[ipg] = volvel;
    }
    _ipg_to_sfc_pressure =
        fz::SafePtr<_FuncRealToCmplx>(_epg_to_sfc_pressure.size());
    for (const auto& [epg, pressure] : _epg_to_sfc_pressure) {
        const _ipg_t ipg = _epg_to_ipg_sfc.at(epg);
        _ipg_to_sfc_pressure[ipg] = pressure;
    }
    _ipg_to_sfc_impedance =
        fz::SafePtr<_FuncRealToCmplx>(_epg_to_sfc_impedance.size());
    for (const auto& [epg, impedance] : _epg_to_sfc_impedance) {
        const _ipg_t ipg = _epg_to_ipg_sfc.at(epg);
        _ipg_to_sfc_impedance[ipg] = impedance;
    }
    _ipg_to_volprop =
        fz::SafePtr<_VolProp>(_epg_to_volprop.size());
    for (const auto& [epg, volprop] : _epg_to_volprop) {
        const _ipg_t ipg = _epg_to_ipg_vol.at(epg);
        _ipg_to_volprop[ipg] = volprop;
    }
    
    // generate the contiguous vector with the ipg of each element
    _ipg_vol_elem = fz::SafePtr<_ipg_t>(_vol_elem_count());
    for (_idx_t e=0; e!=_ipg_vol_elem.size(); ++e) {
        _ipg_vol_elem[e] = _epg_to_ipg_vol.at(_epg_vol_elem[e]);
    }
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::_analize_sparsity()
{
    constexpr std::array<
        std::array<size_t,2>, COMB_REP_SIZE<NODES_IN_3D_ELEM<O>,2>
    > COMB = COMBINATION_REP<NODES_IN_3D_ELEM<O>>;
    _nonzero_row_idx = fz::SafePtr<_idx_t>(_vol_elem_count()*COMB.size());
    _nonzero_col_idx = fz::SafePtr<_idx_t>(_vol_elem_count()*COMB.size());
    auto it_nonzero_row_idx = _nonzero_row_idx.begin();
    auto it_nonzero_col_idx = _nonzero_col_idx.begin();
    for (_idx_t e=0; e!=_vol_elem_count(); ++e) {
        for (size_t i=0; i!=COMB.size(); ++i) {
            *it_nonzero_row_idx = _vol_elem_node_idx[e][COMB[i][0]];
            ++it_nonzero_row_idx;
            *it_nonzero_col_idx = _vol_elem_node_idx[e][COMB[i][1]];
            ++it_nonzero_col_idx;
        }
    }
}

template<ElementOrder O>
Eigen::Matrix<double, NODES_IN_3D_ELEM<O>, 1> shape_func(
    const double&, const double&, const double&
);

template<>
Eigen::Matrix<double, NODES_IN_3D_ELEM<ElementOrder::O1>, 1>
shape_func<ElementOrder::O1>(
    const double& xi0, const double& xi1, const double& xi2
) {
    const double xi3 = 1 - xi0 - xi1 - xi2;
    return Eigen::Matrix<double, NODES_IN_3D_ELEM<ElementOrder::O1>, 1>(
        xi0,
        xi1,
        xi2,
        xi3
    );
}

template<>
Eigen::Matrix<double, NODES_IN_3D_ELEM<ElementOrder::O2>, 1>
shape_func<ElementOrder::O2>(
    const double& xi0, const double& xi1, const double& xi2
) {
    const double xi3 = 1 - xi0 - xi1 - xi2;
    return Eigen::Matrix<double, NODES_IN_3D_ELEM<ElementOrder::O2>, 1>(
        xi0*(2*xi0-1),
        xi1*(2*xi1-1),
        xi2*(2*xi2-1),
        xi3*(2*xi3-1),
        4*xi0*xi1,
        4*xi0*xi2,
        4*xi0*xi3,
        4*xi1*xi2,
        4*xi1*xi3,
        4*xi2*xi3
    );
}

template<ElementOrder O>
Eigen::Matrix<double, DIM, NODES_IN_3D_ELEM<O>> shape_func_gradient(
    const double&, const double&, const double&
);

template<>
Eigen::Matrix<double, DIM, NODES_IN_3D_ELEM<ElementOrder::O1>>
shape_func_gradient<ElementOrder::O1>(
    const double& xi0, const double& xi1, const double& xi2
) {
    return Eigen::Matrix<double, DIM, NODES_IN_3D_ELEM<ElementOrder::O1>> {
        {1, 0, 0, -1},
        {0, 1, 0, -1},
        {0, 0, 1, -1}
    };
}

template<>
Eigen::Matrix<double, DIM, NODES_IN_3D_ELEM<ElementOrder::O2>>
shape_func_gradient<ElementOrder::O2>(
    const double& xi0, const double& xi1, const double& xi2
) {
    const double xi3 = 1 - xi0 - xi1 - xi2;
    return Eigen::Matrix<double, DIM, NODES_IN_3D_ELEM<ElementOrder::O2>> {
        {4*xi0-1, 0, 0, 1-4*xi3, 4*xi1, 4*xi2, 4*(xi3-xi0), 0, -4*xi1, -4*xi2},
        {0, 4*xi1-1, 0, 1-4*xi3, 4*xi0, 0, -4*xi0, 4*xi2, 4*(xi3-xi1), -4*xi2},
        {0, 0, 4*xi2-1, 1-4*xi3, 0, 4*xi0, -4*xi0, 4*xi1, -4*xi1, 4*(xi3-xi2)},
    };
}

template<ElementOrder O>
void SimulationAcFemFreqD3<O>::_solve()
{
    fz::SafePtr<std::complex<double>> ipg_to_density_value(
        _ipg_to_volprop.size()
    );
    fz::SafePtr<std::complex<double>> ipg_to_soundspeed_value(
        _ipg_to_volprop.size()
    );
    for (size_t i=0; i!=_freq_vec.size(); ++i)
    {
        const double freq = _freq_vec[i];
        
        for (_ipg_t ipg=0; ipg!=_ipg_to_volprop.size(); ++ipg) {
            ipg_to_density_value[ipg] = 
                (_ipg_to_volprop[ipg].density)(freq);
            ipg_to_soundspeed_value[ipg] = 
                (_ipg_to_volprop[ipg].soundspeed)(freq);
        }
    }
    ipg_to_density_value.free();
    ipg_to_soundspeed_value.free();
}

template<ElementOrder O> constexpr size_t NGP_STIF = [] {
    if constexpr (O == ElementOrder::O1) { return 1; }
    if constexpr (O == ElementOrder::O2) { return 4; }
}();

template<ElementOrder O> constexpr size_t NGP_MASS = [] {
    if constexpr (O == ElementOrder::O1) { return 4;  }
    if constexpr (O == ElementOrder::O2) { return 15; }
}();

template<size_t N>
constexpr std::array<std::array<double,DIM>,N> GAUSS_POINTS = [] {
    if constexpr (N == 1) {
        return std::array<std::array<double,DIM>,N>({ {1/4, 1/4, 1/4} });
    }
    if constexpr (N == 4) {
        constexpr double a = (5 -   std::sqrt(5)) / 20;
        constexpr double b = (5 + 3*std::sqrt(5)) / 20;
        return std::array<std::array<double,DIM>,N> {{
            {a,a,a}, {b,a,a}, {a,b,a}, {a,a,b}
        }};
    }
    if constexpr (N == 15) {
        constexpr double a = 1/4;
        constexpr double b = ( 7 +   std::sqrt(15)) / 34;
        constexpr double c = ( 7 -   std::sqrt(15)) / 34;
        constexpr double d = (13 - 3*std::sqrt(15)) / 34;
        constexpr double e = (13 + 3*std::sqrt(15)) / 34;
        constexpr double f = ( 5 -   std::sqrt(15)) / 20;
        constexpr double g = ( 5 +   std::sqrt(15)) / 20;
        return std::array<std::array<double,DIM>,N> {{
            {a,a,a}, {b,b,b}, {b,b,d}, {b,d,b}, {d,b,b}, 
            {c,c,c}, {c,c,e}, {c,e,c}, {e,c,c}, {f,f,g},
            {f,g,f}, {g,f,f}, {f,g,g}, {g,f,g}, {g,g,f}
        }};
    }
}();

template<size_t N>
constexpr std::array<double,N> GAUSS_WEIGHTS = [] {
    if constexpr (N == 1) {
        constexpr double a = 1;
        return std::array<double,N>({VOLUME_REF_TET * a});
    }
    if constexpr (N == 4) {
        constexpr double a = VOLUME_REF_TET * 1/4;
        return std::array<double,N>{{a,a,a,a}};
    }
    if constexpr (N == 15) {
        constexpr double a = VOLUME_REF_TET * 8/405;
        constexpr double b = VOLUME_REF_TET * (2665 - 14*std::sqrt(15))/226800;
        constexpr double c = VOLUME_REF_TET * (2665 + 14*std::sqrt(15))/226800;
        constexpr double d = VOLUME_REF_TET * 5/567;
        return std::array<double,N>{{a,b,b,b,b,c,c,c,c,d,d,d,d,d,d}};
    }
}();

template <ElementOrder O>
ResultAcFemFreqD3 SimulationAcFemFreqD3<O>::run()
{
    _check_if_it_can_run();
    _define_freq_vector();
    _organize_physical_group_data();
    _analize_sparsity();
    
    // matrices
    constexpr std::array<
        std::array<size_t,2>, COMB_REP_SIZE<NODES_IN_3D_ELEM<O>,2>
    > COMB = COMBINATION_REP<NODES_IN_3D_ELEM<O>>;
    _btb_detj_w_vals = fz::SafePtr<_idx_t>(_vol_elem_count()*COMB.size());
    _btb_detj_w_vals.fill(0);
    auto it_btb_detj_w_vals = _btb_detj_w_vals.begin();

    for (_idx_t e=0; e!=_vol_elem_count(); ++e)
    {
        // coordinates matrix
        Eigen::Matrix<double,NODES_IN_3D_ELEM<O>,DIM> coords_matrix;
        for (size_t i=0; i!=NODES_IN_3D_ELEM<O>; ++i) {
            const _idx_t node_idx = _vol_elem_node_idx[e][i];
            for (size_t j=0; j!=DIM; ++j) {
                coords_matrix(i,j) = _node_coords[node_idx][j];
            }
        }

        // stiffness matrices
        constexpr std::array<std::array<double,DIM>,NGP_STIF<O>>
            gauss_points = GAUSS_POINTS<NGP_STIF<O>>;
        it_btb_detj_w_vals += COMB_REP_SIZE<NODES_IN_3D_ELEM<O>,2>;
        for (_idx_t p=0; p!=NGP_STIF<O>; ++p)
        {
            Eigen::Matrix<double,DIM,NODES_IN_3D_ELEM<O>> nabla_n =
                shape_func_gradient<O>(
                    gauss_points[p][0], gauss_points[p][1], gauss_points[p][2]
                ); // todo: try putting constexpr here

            const Eigen::Matrix<double,DIM,DIM> jac_matrix =
                nabla_n * coords_matrix;

            const Eigen::Matrix<double,DIM,DIM> inv_jac = jac_matrix.inverse();
            
            const Eigen::Matrix<double,DIM,NODES_IN_3D_ELEM<O>> b_matrix =
                inv_jac * nabla_n;
            
            const Eigen::Matrix<double,NODES_IN_3D_ELEM<O>,NODES_IN_3D_ELEM<O>>
                btb = b_matrix.transpose() * b_matrix;
            
            const double det_jac = jac_matrix.determinant();
            
            const Eigen::Matrix<double,NODES_IN_3D_ELEM<O>,NODES_IN_3D_ELEM<O>>
                btb_detj_w = btb * det_jac * GAUSS_WEIGHTS<NGP_STIF<O>>[p];
            
            for (size_t c=0; c!=COMB_REP_SIZE<NODES_IN_3D_ELEM<O>,2>; ++c) {
                *it_btb_detj_w_vals += btb_detj_w(COMB[c][0], COMB[c][1]);
                ++it_btb_detj_w_vals;
            }
            it_btb_detj_w_vals -= COMB_REP_SIZE<NODES_IN_3D_ELEM<O>,2>;
        }
    }
    
    

    _solve();
    return ResultAcFemFreqD3();
}

void trim_right_whitespace(std::string_view& sv) {
    constexpr std::string_view WHITE_SPACE = " \t\n\r\f\v";
    const size_t end = sv.find_last_not_of(WHITE_SPACE);
    sv = (end == std::string_view::npos) ? "" : sv.substr(0, end+1);
}

template<typename T>
T parse(std::string_view str) {
    trim_right_whitespace(str);
    T value;
    auto result = std::from_chars(str.data(), str.data() + str.size(), value);
    if (result.ec != std::errc{} || result.ptr != str.data()+str.size()) {
        throw std::invalid_argument("invalid number format");
    }
    return value;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::_load_bdf(const char* const path_to_mesh)
{
    constexpr size_t MAX_BDF_CHARACTERS_PER_LINE = 80;
    std::ifstream file(path_to_mesh);
    std::string line;
    line.reserve(MAX_BDF_CHARACTERS_PER_LINE);
    
    if (!file.is_open()) {
        log::error("Could not open file: {}", path_to_mesh);
    }

    // first pass: count lines by type
    size_t node_count = 0;
    size_t sfc_elem_count = 0;
    size_t vol_elem_count = 0;
    while (std::getline(file, line)) {
        if      (line.starts_with("GRID"))    { ++node_count; }
        else if (line.starts_with("CTRIA3"))  { ++sfc_elem_count; }
        else if (line.starts_with("CTETRA"))  { ++vol_elem_count; }
        else if (line.starts_with("ENDDATA")) { break; }
    }

    // second pass: parse data
    _node_coords = fz::SafePtr<std::array<double,3>>(node_count);
    _sfc_elem_node_idx =
        fz::SafePtr<std::array<size_t,NODES_IN_2D_ELEM<O>>>(sfc_elem_count);
    _vol_elem_node_idx =
        fz::SafePtr<std::array<size_t,NODES_IN_3D_ELEM<O>>>(vol_elem_count);
    _epg_sfc_elem = fz::SafePtr<size_t>(sfc_elem_count);
    _epg_vol_elem = fz::SafePtr<size_t>(vol_elem_count);
    auto it_node_coords = _node_coords.begin();
    auto it_2d_elem_vtx_idx = _sfc_elem_node_idx.begin();
    auto it_3d_elem_vtx_idx = _vol_elem_node_idx.begin();
    auto it_2d_elem_epg = _epg_sfc_elem.begin();
    auto it_3d_elem_epg = _epg_vol_elem.begin();
    file.clear();
    file.seekg(0, std::ios::beg);
    while (std::getline(file, line)) {
        if (line.starts_with("GRID")) {
            *it_node_coords = {
                parse<double>(line.substr(24,8)),
                parse<double>(line.substr(32,8)),
                parse<double>(line.substr(40,8))
            };
            ++it_node_coords;
        }
        else if (line.starts_with("CTRIA3")) {
            const _epg_t epg = parse<size_t>(line.substr(16,8));
            _existing_epg_sfc.insert(epg);
            *it_2d_elem_epg = epg;
            ++it_2d_elem_epg;
            *it_2d_elem_vtx_idx = { // minus one for zero base indexing
                parse<size_t>(line.substr(24,8)) - 1,
                parse<size_t>(line.substr(32,8)) - 1,
                parse<size_t>(line.substr(40,8)) - 1
            };
            ++it_2d_elem_vtx_idx;
        }
        else if (line.starts_with("CTETRA")) {
            const _epg_t epg = parse<size_t>(line.substr(16,8));
            _existing_epg_vol.insert(epg);
            *it_3d_elem_epg = epg;
            ++it_3d_elem_epg;
            *it_3d_elem_vtx_idx = { // minus one for zero base indexing
                parse<size_t>(line.substr(24,8)) - 1,
                parse<size_t>(line.substr(32,8)) - 1,
                parse<size_t>(line.substr(40,8)) - 1,
                parse<size_t>(line.substr(48,8)) - 1
            };
            ++it_3d_elem_vtx_idx;
        }
        else if (line.starts_with("ENDDATA")) {
            break;
        }
    }
    file.close();
}

template<>
void SimulationAcFemFreqD3<ElementOrder::O1>::_generate_extra_nodes() {
    // nothing needs to be done in this case (in order 1)
}

template<>
void SimulationAcFemFreqD3<ElementOrder::O2>::_generate_extra_nodes()
{
    constexpr std::array<
        std::array<size_t,2>,EXTRA_NODES_IN_3D_ELEM<ElementOrder::O2>
    > VTX_PAIRS_3D = {{ {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3} }};

    constexpr std::array<
        std::array<size_t,2>,EXTRA_NODES_IN_2D_ELEM<ElementOrder::O2>
    > VTX_PAIRS_2D = {{ {0,1}, {0,2}, {1,2} }};

    std::unordered_map<std::tuple<size_t,size_t>,size_t> idxs_extra_nodes;
    
    // first pass: count extra nodes and save idx tuples
    fz::SafePtr<std::array<bool,EXTRA_NODES_IN_3D_ELEM<ElementOrder::O2>>>
        is_extra_node(_vol_elem_count());
    size_t count = _node_count();
    for (size_t e=0; e!=_vol_elem_count(); ++e) {
        for (size_t i=0; i!=VTX_PAIRS_3D.size(); ++i)
        {
            const std::tuple<size_t,size_t> tup = make_ascending_tuple(
                _vol_elem_node_idx[e][VTX_PAIRS_3D[i][0]],
                _vol_elem_node_idx[e][VTX_PAIRS_3D[i][1]]
            );
            if (!idxs_extra_nodes.contains(tup)) {
                is_extra_node[e][i] = true;
                _vol_elem_node_idx[e][NODES_IN_3D_ELEM<ElementOrder::O1> + i] =
                    count;
                idxs_extra_nodes.insert({tup, count});
                ++count;
            } else {
                is_extra_node[e][i] = false;          
                _vol_elem_node_idx[e][NODES_IN_3D_ELEM<ElementOrder::O1> + i] =
                    idxs_extra_nodes.at(tup);
            }     
        }
    }

    // TODO: grow() here
    auto temp = std::move(_node_coords);
    _node_coords = fz::SafePtr<std::array<double,3>>(count);
    std::copy(temp.begin(), temp.end(), _node_coords.begin());
    temp.free();
    
    // second pass: create the extra nodes and assign to 3D elements
    for (size_t e=0; e!=_vol_elem_count(); ++e) {
        for (size_t i=0; i!=VTX_PAIRS_3D.size(); ++i)
        {   
            if (!is_extra_node[e][i]) { continue; }

            const std::tuple<size_t,size_t> tup = make_ascending_tuple(
                _vol_elem_node_idx[e][VTX_PAIRS_3D[i][0]],
                _vol_elem_node_idx[e][VTX_PAIRS_3D[i][1]]
            );
            const double x = mean(
                _node_coords[std::get<0>(tup)][0],
                _node_coords[std::get<1>(tup)][0]
            );
            const double y = mean(
                _node_coords[std::get<0>(tup)][1],
                _node_coords[std::get<1>(tup)][1]
            );
            const double z = mean(
                _node_coords[std::get<0>(tup)][2],
                _node_coords[std::get<1>(tup)][2]
            );
            const size_t idx_extra_node = idxs_extra_nodes.at(tup);
            _node_coords[idx_extra_node] = {x, y, z}; // add extra node
        }
    }
    is_extra_node.free();

    // third pass: assign nodes to 2D elements
    for (size_t e=0; e!=_sfc_elem_count(); ++e) {
        for (size_t i=0; i!=VTX_PAIRS_2D.size(); ++i)
        {
            // create a tuple of the indices in ascending order
            const std::tuple<size_t,size_t> tup = make_ascending_tuple(
                _sfc_elem_node_idx[e][VTX_PAIRS_2D[i][0]],
                _sfc_elem_node_idx[e][VTX_PAIRS_2D[i][1]]
            );
            _sfc_elem_node_idx[e][NODES_IN_2D_ELEM<ElementOrder::O1> + i] =
                idxs_extra_nodes.at(tup);
        }
    }
}

} // namespace numav
