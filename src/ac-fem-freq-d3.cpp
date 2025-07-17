// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav.hpp"
#include "logger.hpp"

inline numav::Logger _logger;

#include <tuple>

#include <charconv>
#include <fstream>
#undef NDEBUG
#include <cassert>

#include "Eigen/Dense"

#include "typedefs.hpp"

using ResultAcFemFreqD3 = typename numav::Result<
    numav::Phenomenon::ACOUSTIC,
    numav::NumericalMethod::FEM,
    numav::Domain::FREQUENCY,
    numav::Dimension::D3
>;

ResultAcFemFreqD3::Result() = default;
ResultAcFemFreqD3::~Result() = default;

template<numav::ElementOrder O>
using SimulationAcFemFreqD3 = typename numav::Simulation<
    numav::Phenomenon::ACOUSTIC,
    numav::NumericalMethod::FEM,
    numav::Domain::FREQUENCY,
    numav::Dimension::D3,
    O
>;

template<numav::ElementOrder O>
SimulationAcFemFreqD3<O>::Simulation() {
    _is_mesh_defined = false;
    _is_freq_range_defined = false;
    _is_any_source_defined = false;
}

template<numav::ElementOrder O>
SimulationAcFemFreqD3<O>::~Simulation() {
    _node_coords.free();
    _sfc_elem_node_idx.free();
    _vol_elem_node_idx.free();
    _epg_sfc_elem.free();
    _epg_vol_elem.free();
    _freq_vec.free();
}

fz::SafePtr<double> linspace(
    const double& start,
    const double& finish,
    const size_t& num_points
) {
    assert(num_points!=0 && num_points!=1);
    fz::SafePtr<double> result(num_points);
    double step = (finish - start) / static_cast<double>(num_points - 1);
    result.front() = start;
    for (double* it=result.begin()+1; it!=result.end()-1; ++it) {
        *it = *(it-1) + step;
    }
    result.back() = finish;
    return result;
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::_check_if_mesh_is_defined() {
    if (!_is_mesh_defined){
        _logger.error("Mesh not defined. Call load_mesh to do so.");
    }
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::set_freq_range(
    const double& freq_min,
    const double& freq_max
) {
    _freq_min = freq_min;
    _freq_max = freq_max;
    _is_freq_range_defined = true;
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::load_mesh(
    const char* const path_to_mesh
) {
    if (_is_mesh_defined) {
        _logger.error("Mesh is already defined.");
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
        _logger.error("Unrecognized file format: \"{0}\"", format);
    }
    _is_mesh_defined = true;
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::add_volume_material(
    const _epg_t& epg,
    const _FuncRealToCmplx& density,
    const _FuncRealToCmplx& soundspeed
) {
    _check_if_mesh_is_defined();
    if (!_existing_epg_vol.contains(epg)) {
        // TODO: error
    }
    if (_epg_to_volprop.contains(epg)) {
        // TODO: error
    }
    _epg_to_volprop.insert({epg, {density, soundspeed}});
    const _ipg_t ipg = _epg_to_ipg_vol.size();
    _epg_to_ipg_vol.insert({epg, ipg});
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::add_source(
    const TypeOfSource& type_of_source,
    const std::array<double,3>& point_coordinates,
    const PhysicalQuantity& physical_quantity_type,
    const _FuncRealToCmplx& physical_quantity_value
) {
    _check_if_mesh_is_defined();
    if (type_of_source != numav::TypeOfSource::POINT) {
        // TODO: error
    }
    const _idx_t closest_point_idx = _get_closest_point(point_coordinates);
    
    if (physical_quantity_type == numav::PhysicalQuantity::VOLUME_VELOCITY) {
        _point_volvel.push_back(
            std::make_tuple(closest_point_idx, physical_quantity_value)
        );
    }
    else if (physical_quantity_type == numav::PhysicalQuantity::PRESSURE) {
        _point_volvel.push_back(
            std::make_tuple(closest_point_idx, physical_quantity_value)
        ); 
    }
    else {
        // TODO: error
    }
    _is_any_source_defined = true;
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::add_source(
    const TypeOfSource& type_of_source,
    const _epg_t& epg,
    const PhysicalQuantity& physical_quantity_type,
    const _FuncRealToCmplx& physical_quantity_value
) {
    _check_if_mesh_is_defined();
    if (type_of_source != numav::TypeOfSource::SURFACE) {
        // TODO: error
    }
    if (!_existing_epg_sfc.contains(epg)) {
        // TODO: error
    }
    if (_epg_to_sfc_pressure.contains(epg) ||
        _epg_to_sfc_volvel.contains(epg) ||
        _epg_to_sfc_impedance.contains(epg)
    ) {
        // TODO: error
    }
    if (physical_quantity_type == numav::PhysicalQuantity::PRESSURE) {
        _epg_to_sfc_pressure.insert({epg, physical_quantity_value});
        const _ipg_t ipg = _epg_to_ipg_sfc.size();
        _epg_to_ipg_sfc.insert({epg, ipg});
        _is_any_source_defined = true;
        return;
    }
    if (physical_quantity_type == numav::PhysicalQuantity::VOLUME_VELOCITY) {
        _epg_to_sfc_volvel.insert({epg, physical_quantity_value});
        const _ipg_t ipg = _epg_to_ipg_sfc.size();
        _epg_to_ipg_sfc.insert({epg, ipg});
        _is_any_source_defined = true;
        return;
    }
    // TODO: error
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::add_surface_specific_acoustic_impedance(
    const _epg_t& epg,
    const _FuncRealToCmplx& impedance
) {
    _check_if_mesh_is_defined();
    if (!_existing_epg_sfc.contains(epg)) {
        // TODO: error
    }
    if (_epg_to_sfc_pressure.contains(epg) ||
        _epg_to_sfc_volvel.contains(epg) ||
        _epg_to_sfc_impedance.contains(epg)
    ) {
        // TODO: error
    }
    _epg_to_sfc_impedance.insert({epg, impedance});
    const _ipg_t ipg = _epg_to_ipg_sfc.size();
    _epg_to_ipg_sfc.insert({epg, ipg});
}

template <numav::ElementOrder O>
_idx_t SimulationAcFemFreqD3<O>::_get_closest_point(
    const std::array<double,3>&
) {
    return 0; // TODO
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::_check_if_it_can_run() {
    _check_if_mesh_is_defined();
    if (!_is_any_source_defined){
        // TODO error
    }
    if (!_is_freq_range_defined){
        // TODO error
    }
    for (auto& epg : _existing_epg_vol) {
        if (!_epg_to_volprop.contains(epg)) {
            // TODO: error
        }
    }
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::_define_freq_vector() {
    // TODO: decide number here
    // TODO: make it not linear
    _freq_vec = linspace(_freq_min, _freq_max, 8192);
}

template <numav::ElementOrder O>
ResultAcFemFreqD3 SimulationAcFemFreqD3<O>::run()
{
    _check_if_it_can_run();
    _define_freq_vector();

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
    
    _ipg_vol_elem = fz::SafePtr<_ipg_t>(_vol_elem_count); // todo

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

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::_load_bdf(const char* const path_to_mesh)
{
    constexpr size_t MAX_BDF_CHARACTERS_PER_LINE = 80;
    std::ifstream file(path_to_mesh);
    std::string line;
    line.reserve(MAX_BDF_CHARACTERS_PER_LINE);
    
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << path_to_mesh << "\n";
        return;
    }

    // first pass: count lines by type
    _node_count = 0;
    _sfc_elem_count = 0;
    _vol_elem_count = 0;
    while (std::getline(file, line)) {
        if      (line.starts_with("GRID"))   { ++_node_count;    }
        else if (line.starts_with("CTRIA3")) { ++_sfc_elem_count; }
        else if (line.starts_with("CTETRA")) { ++_vol_elem_count; }
    }

    // second pass: parse data
    _node_coords = fz::SafePtr<std::array<double,3>>(_node_count);
    _sfc_elem_node_idx = fz::SafePtr<std::array<size_t,NODES_IN_2D_ELEM<O>>>(
        _sfc_elem_count
    );
    _vol_elem_node_idx = fz::SafePtr<std::array<size_t,NODES_IN_3D_ELEM<O>>>(
        _vol_elem_count
    );
    _epg_sfc_elem = fz::SafePtr<size_t>(_sfc_elem_count);
    _epg_vol_elem = fz::SafePtr<size_t>(_vol_elem_count);
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
    }
    file.close();
    _generate_extra_nodes(); // call is based on the element order
}

template<>
void SimulationAcFemFreqD3<numav::ElementOrder::O1>::_generate_extra_nodes() {
    // nothing needs to be done for this case
}

// define a hash function for std::tuple<size_t,size_t>
namespace std {
    template<>
    struct hash<std::tuple<size_t,size_t>> {
        size_t operator()(const std::tuple<size_t,size_t>& key) const {
            return (std::get<0>(key) << 4*sizeof(size_t)) + std::get<1>(key);
        }
    };
}

template<typename... T>
auto mean(const T&... args) {
    return (args + ...) / (sizeof...(args));
}

template<typename T>
std::tuple<T,T> make_ascending_tuple(const T& a, const T& b) {
    return a<b ? std::make_tuple(a,b) : std::make_tuple(b,a);
}

template<>
void SimulationAcFemFreqD3<numav::ElementOrder::O2>::_generate_extra_nodes()
{
    constexpr std::array<
        std::array<size_t,2>,EXTRA_NODES_IN_3D_ELEM<ElementOrder::O2>
    > VTX_PAIRS_3D = {{ {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3} }};

    constexpr std::array<
        std::array<size_t,2>,EXTRA_NODES_IN_2D_ELEM<ElementOrder::O2>
    > VTX_PAIRS_2D = {{ {0,1}, {0,2}, {1,2} }};

    std::unordered_map<std::tuple<size_t,size_t>,size_t> idxs_extra_nodes;
    
    // first pass: count extra nodes and save idx tuples
    fz::SafePtr<std::array<
        bool, EXTRA_NODES_IN_3D_ELEM<ElementOrder::O2>
    >> is_extra_node(_vol_elem_count);
    for (size_t e=0; e!=_vol_elem_count; ++e) {
        for (size_t i=0; i!=VTX_PAIRS_3D.size(); ++i)
        {
            const std::tuple<size_t,size_t> tup = make_ascending_tuple(
                _vol_elem_node_idx[e][VTX_PAIRS_3D[i][0]],
                _vol_elem_node_idx[e][VTX_PAIRS_3D[i][1]]
            );
            if (!idxs_extra_nodes.contains(tup)) {
                is_extra_node[e][i] = true;
                _vol_elem_node_idx[e][NODES_IN_3D_ELEM<ElementOrder::O1> + i] =
                    _node_count;
                idxs_extra_nodes.insert({tup, _node_count});
                ++_node_count;
            } else {
                is_extra_node[e][i] = false;          
                _vol_elem_node_idx[e][NODES_IN_3D_ELEM<ElementOrder::O1> + i] =
                    idxs_extra_nodes.find(tup)->second;
            }     
        }
    }

    // TODO: grow() here
    auto temp = std::move(_node_coords);
    _node_coords = fz::SafePtr<std::array<double,3>>(_node_count);
    std::copy(temp.begin(), temp.end(), _node_coords.begin());
    temp.free();
    
    // second pass: create the extra nodes and assign to 3D elements
    for (size_t e=0; e!=_vol_elem_count; ++e) {
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
            const size_t idx_extra_node = idxs_extra_nodes.find(tup)->second;
            _node_coords[idx_extra_node] = {x, y, z}; // add extra node
        }
    }
    is_extra_node.free();

    // third pass: assign nodes to 2D elements
    for (size_t e=0; e!=_sfc_elem_count; ++e) {
        for (size_t i=0; i!=VTX_PAIRS_2D.size(); ++i)
        {
            // create a tuple of the indices in ascending order
            const std::tuple<size_t,size_t> tup = make_ascending_tuple(
                _sfc_elem_node_idx[e][VTX_PAIRS_2D[i][0]],
                _sfc_elem_node_idx[e][VTX_PAIRS_2D[i][1]]
            );
            _sfc_elem_node_idx[e][NODES_IN_2D_ELEM<ElementOrder::O1> + i] =
                idxs_extra_nodes.find(tup)->second;
        }
    }
}
