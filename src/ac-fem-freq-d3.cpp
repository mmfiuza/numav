// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav.hpp"

#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <charconv>
#include <fstream>
#undef NDEBUG
#include <cassert>

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
    _is_freq_defined = false;
    _is_any_source_defined = false;
}

template<numav::ElementOrder O>
SimulationAcFemFreqD3<O>::~Simulation() {
    _freq_vector.free();
}

SafePtr<double> linspace(
    const double& start, const double& finish, const size_t& num_points
) {
    assert(num_points!=0 && num_points!=1);
    SafePtr<double> result(num_points);
    double step = (finish - start) / static_cast<double>(num_points - 1);
    result.front() = start;
    for (double* it=result.begin()+1; it!=result.end()-1; ++it) {
        *it = *(it-1) + step;
    }
    result.back() = finish;
    return result;
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::set_freq_limits(
    const double& freq_min, const double& freq_max
) {
    _freq_min = freq_min;
    _freq_max = freq_max;
    // TODO: decide number here
    _freq_vector = linspace(_freq_min, _freq_max, 8192);
    _is_freq_defined = true;
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::load_mesh(const char* const path_to_mesh) {
    _mesh.load_bdf(path_to_mesh);
    _is_mesh_defined = true;
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::add_volume_material(
    const uint64_t& id, const double& rho, const double& c
) {
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::add_source(
    const TypeOfSource& type_of_source,
    const std::array<double,3>& point_coordinates,
    const PhysicalQuantity& physical_quantity,
    const std::function<std::complex<double>(double)>& physical_quantity_value
) {
    _is_any_source_defined = true;
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::add_source(
    const TypeOfSource& type_of_source,
    const uint64_t& surface_id,
    const PhysicalQuantity& physical_quantity,
    const std::function<std::complex<double>(double)>& physical_quantity_value
) {
    _is_any_source_defined = true;
}

template <numav::ElementOrder O>
void SimulationAcFemFreqD3<O>::add_surface_specific_acoustic_impedance(
    const uint64_t& surface_id,
    const std::function<std::complex<double>(double)>& impedance
) {
}

template <numav::ElementOrder O>
ResultAcFemFreqD3 SimulationAcFemFreqD3<O>::run() {
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
void SimulationAcFemFreqD3<O>::Mesh::load_bdf(const char* const path_to_mesh)
{
    constexpr size_t MAX_BDF_CHARACTERS_PER_LINE = 80;
    std::ifstream file(path_to_mesh);
    std::string line;
    line.reserve(MAX_BDF_CHARACTERS_PER_LINE);
    
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << path_to_mesh << "\n";
        return;
    }

    _node_count = 0;
    _2d_elem_count = 0;
    _3d_elem_count = 0;
    _next_2d_elem_tag = 0;
    _next_3d_elem_tag = 0;
    
    // first pass: count lines by type and fill the file_id_to_tag maps
    while (std::getline(file, line)) {
        if (line.starts_with("GRID")) {
            ++_node_count;
        }
        else if (line.starts_with("CTRIA3")) {
            size_t elem_file_id = parse<size_t>(line.substr(16,8));
            if (!_file_id_to_tag_2d.contains(elem_file_id)) {
                size_t elem_tag = _next_2d_elem_tag;
                ++_next_2d_elem_tag;
                _file_id_to_tag_2d.insert({elem_file_id, elem_tag});
            }
            ++_2d_elem_count;
        }
        else if (line.starts_with("CTETRA")) {
            size_t elem_file_id = parse<size_t>(line.substr(16,8));
            if (!_file_id_to_tag_3d.contains(elem_file_id)) {
                size_t elem_tag = _next_3d_elem_tag;
                ++_next_3d_elem_tag;
                _file_id_to_tag_3d.insert({elem_file_id, elem_tag});
            }
            ++_3d_elem_count;
        }
    }

    // second pass: parse data
    _node_coords = SafePtr<std::array<double,DIM_COUNT<Dimension::D3>>>(
        _node_count
    );
    _2d_elem_vtx_idx = SafePtr<std::array<size_t,NODES_IN_2D_ELEM<O>>>(
        _2d_elem_count
    );
    _3d_elem_vtx_idx = SafePtr<std::array<size_t,NODES_IN_3D_ELEM<O>>>(
        _3d_elem_count
    );
    _2d_elem_tag = SafePtr<size_t>(_2d_elem_count);
    _3d_elem_tag = SafePtr<size_t>(_3d_elem_count);
    auto it_node_coords = _node_coords.begin();
    auto it_2d_elem_vtx_idx = _2d_elem_vtx_idx.begin();
    auto it_3d_elem_vtx_idx = _3d_elem_vtx_idx.begin();
    auto it_2d_elem_tag = _2d_elem_tag.begin();
    auto it_3d_elem_tag = _3d_elem_tag.begin();
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
            *it_2d_elem_tag = _file_id_to_tag_2d.find(
                parse<size_t>(line.substr(16,8))
            )->second;
            ++it_2d_elem_tag;
            *it_2d_elem_vtx_idx = {
                parse<size_t>(line.substr(24,8)),
                parse<size_t>(line.substr(32,8)),
                parse<size_t>(line.substr(40,8))
            };
            ++it_2d_elem_vtx_idx;
        }
        else if (line.starts_with("CTETRA")) {
            *it_3d_elem_tag = _file_id_to_tag_3d.find(
                parse<size_t>(line.substr(16,8))
            )->second;
            ++it_3d_elem_tag;
            *it_3d_elem_vtx_idx = {
                parse<size_t>(line.substr(24,8)),
                parse<size_t>(line.substr(32,8)),
                parse<size_t>(line.substr(40,8)),
                parse<size_t>(line.substr(48,8))
            };
            ++it_3d_elem_vtx_idx;
        }
    }
    file.close();
    _generate_extra_nodes(); // call is based on the element order
}

template<>
void SimulationAcFemFreqD3<numav::ElementOrder::O1>::
Mesh::_generate_extra_nodes() {
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
void SimulationAcFemFreqD3<numav::ElementOrder::O2>::
Mesh::_generate_extra_nodes()
{
    constexpr std::array<
        std::array<size_t,2>,EXTRA_NODES_IN_3D_ELEM<ElementOrder::O2>
    > VTX_PAIRS_3D = {{ {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3} }};

    constexpr std::array<
        std::array<size_t,2>,EXTRA_NODES_IN_2D_ELEM<ElementOrder::O2>
    > VTX_PAIRS_2D = {{ {0,1}, {0,2}, {1,2} }};

    std::unordered_map<std::tuple<size_t,size_t>,size_t> idxs_extra_nodes;
    
    // first pass: count extra nodes and save idx tuples
    SafePtr<std::array<
        bool, EXTRA_NODES_IN_3D_ELEM<ElementOrder::O2>
    >> is_extra_node(_3d_elem_count);
    for (size_t e=0; e!=_3d_elem_count; ++e) {
        for (size_t i=0; i!=VTX_PAIRS_3D.size(); ++i)
        {
            const std::tuple<size_t,size_t> tup = make_ascending_tuple(
                _3d_elem_vtx_idx[e][VTX_PAIRS_3D[i][0]],
                _3d_elem_vtx_idx[e][VTX_PAIRS_3D[i][1]]
            );
            if (!idxs_extra_nodes.contains(tup)) {
                is_extra_node[e][i] = true;
                _3d_elem_vtx_idx[e][NODES_IN_3D_ELEM<ElementOrder::O1> + i] =
                    _node_count;
                idxs_extra_nodes.insert({tup, _node_count});
                ++_node_count;
            } else {
                is_extra_node[e][i] = false;          
                _3d_elem_vtx_idx[e][NODES_IN_3D_ELEM<ElementOrder::O1> + i] =
                    idxs_extra_nodes.find(tup)->second;
            }     
        }
    }

    // TODO: grow() here
    auto temp = std::move(_node_coords);
    _node_coords =
        SafePtr<std::array<double,DIM_COUNT<Dimension::D3>>>(_node_count);
    std::copy(temp.begin(), temp.end(), _node_coords.begin());
    temp.free();
    
    // second pass: create the extra nodes and assign to 3D elements
    for (size_t e=0; e!=_3d_elem_count; ++e) {
        for (size_t i=0; i!=VTX_PAIRS_3D.size(); ++i)
        {   
            if (!is_extra_node[e][i]) { continue; }

            const std::tuple<size_t,size_t> tup = make_ascending_tuple(
                _3d_elem_vtx_idx[e][VTX_PAIRS_3D[i][0]],
                _3d_elem_vtx_idx[e][VTX_PAIRS_3D[i][1]]
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
    for (size_t e=0; e!=_2d_elem_count; ++e) {
        for (size_t i=0; i!=VTX_PAIRS_2D.size(); ++i)
        {
            // create a tuple of the indices in ascending order
            const std::tuple<size_t,size_t> tup = make_ascending_tuple(
                _2d_elem_vtx_idx[e][VTX_PAIRS_2D[i][0]],
                _2d_elem_vtx_idx[e][VTX_PAIRS_2D[i][1]]
            );
            _2d_elem_vtx_idx[e][NODES_IN_2D_ELEM<ElementOrder::O1> + i] =
                idxs_extra_nodes.find(tup)->second;
        }
    }
}

template<numav::ElementOrder O>
SimulationAcFemFreqD3<O>::Simulation::Mesh::Mesh() {
};

template<numav::ElementOrder O>
SimulationAcFemFreqD3<O>::Simulation::Mesh::~Mesh() {
    _node_coords.free();
    _2d_elem_vtx_idx.free();
    _3d_elem_vtx_idx.free();
    _2d_elem_tag.free();
    _3d_elem_tag.free();
};
