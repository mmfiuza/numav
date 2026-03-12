// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "modules/ac-fem-freq-d3/impl.hpp"

#include <tuple>
#include <fstream>
#include <limits>

#include "common/exception.hpp"
#include "common/log.hpp"
#include "common/hash-functions.hpp"
#include "common/maths.hpp"
#include "common/utils.hpp"

namespace numav {

template <ElementOrder O>
size_t SimulationAcFemFreqD3<O>::Impl::_get_closest_point(
    const std::array<double,3>& point_coords
) {
    double minimum_distance_squared = std::numeric_limits<double>::max();
    size_t ni_closest = std::numeric_limits<size_t>::max();
    for (size_t ni=0; ni!=_ni_count; ++ni) {
        double distance_squared = 
            std::pow(_ni_to_coords[ni][0] - point_coords[0], 2) +
            std::pow(_ni_to_coords[ni][1] - point_coords[1], 2) +
            std::pow(_ni_to_coords[ni][2] - point_coords[2], 2)
        ;
        if (distance_squared < minimum_distance_squared) {
            minimum_distance_squared = distance_squared;
            ni_closest = ni;
        }
    }
    return ni_closest;
}

template <ElementOrder O>
void SimulationAcFemFreqD3<O>::Impl::_load_bdf(const char* const path_to_mesh)
{
    constexpr size_t MAX_BDF_CHARACTERS_PER_LINE = 80;
    std::ifstream file(path_to_mesh);
    std::string line;
    line.reserve(MAX_BDF_CHARACTERS_PER_LINE);
    
    if (!file.is_open()) {
        error("Could not open file: {}", path_to_mesh);
    }

    // first pass: count lines by type
    _ni_count = 0;
    _sei_count = 0;
    _vei_count = 0;
    while (std::getline(file, line)) {
        if      (line.starts_with("GRID"))    { ++_ni_count; }
        else if (line.starts_with("CTRIA3"))  { ++_sei_count; }
        else if (line.starts_with("CTETRA"))  { ++_vei_count; }
        else if (line.starts_with("ENDDATA")) { break; }
    }

    // second pass: parse data
    _ni_to_coords = fz::SafePtr<std::array<double,3>>(_ni_count);
    _sei_to_ni =
        fz::SafePtr<std::array<size_t,NODES_IN_SFC_ELEM<O>>>(_sei_count);
    _vei_to_ni =
        fz::SafePtr<std::array<size_t,NODES_IN_VOL_ELEM<O>>>(_vei_count);
    _sei_to_espg = fz::SafePtr<size_t>(_sei_count);
    _vei_to_evpg = fz::SafePtr<size_t>(_vei_count);
    auto it_ni_to_coords = _ni_to_coords.begin();
    auto it_sfc_elem_node_idx = _sei_to_ni.begin();
    auto it_vol_elem_node_idx = _vei_to_ni.begin();
    auto it_elem_idx_to_espg = _sei_to_espg.begin();
    auto it_elem_idx_to_evpg = _vei_to_evpg.begin();
    file.clear();
    file.seekg(0, std::ios::beg);
    while (std::getline(file, line)) {
        if (line.starts_with("GRID")) {
            *it_ni_to_coords = {
                parse<double>(line.substr(24,8)),
                parse<double>(line.substr(32,8)),
                parse<double>(line.substr(40,8))
            };
            ++it_ni_to_coords;
        }
        else if (line.starts_with("CTRIA3")) {
            const size_t espg = parse<size_t>(line.substr(16,8));
            _existing_espg.insert(espg);
            *it_elem_idx_to_espg = espg;
            ++it_elem_idx_to_espg;
            *it_sfc_elem_node_idx = { // minus one for zero base indexing
                parse<size_t>(line.substr(24,8)) - 1,
                parse<size_t>(line.substr(32,8)) - 1,
                parse<size_t>(line.substr(40,8)) - 1
            };
            ++it_sfc_elem_node_idx;
        }
        else if (line.starts_with("CTETRA")) {
            const size_t evpg = parse<size_t>(line.substr(16,8));
            _existing_evpg.insert(evpg);
            *it_elem_idx_to_evpg = evpg;
            ++it_elem_idx_to_evpg;
            *it_vol_elem_node_idx = { // minus one for zero base indexing
                parse<size_t>(line.substr(24,8)) - 1,
                parse<size_t>(line.substr(32,8)) - 1,
                parse<size_t>(line.substr(40,8)) - 1,
                parse<size_t>(line.substr(48,8)) - 1
            };
            ++it_vol_elem_node_idx;
        }
        else if (line.starts_with("ENDDATA")) {
            break;
        }
    }
    file.close();
}

template<>
void SimulationAcFemFreqD3<ElementOrder::O1>::Impl::_generate_extra_nodes() {
    // nothing needs to be done in this case (in order 1)
}

template<>
void SimulationAcFemFreqD3<ElementOrder::O2>::Impl::_generate_extra_nodes()
{
    constexpr std::array<
        std::array<size_t,2>,EXTRA_NODES_IN_VOL_ELEM<ElementOrder::O2>
    > VTX_PAIRS_VOL = {{ {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3} }};

    constexpr std::array<
        std::array<size_t,2>,EXTRA_NODES_IN_SFC_ELEM<ElementOrder::O2>
    > VTX_PAIRS_SFC = {{ {0,1}, {0,2}, {1,2} }};

    std::unordered_map<std::tuple<size_t,size_t>,size_t> idxs_extra_nodes;
    
    // first pass: count extra nodes and save idx tuples
    fz::SafePtr<std::array<bool,EXTRA_NODES_IN_VOL_ELEM<ElementOrder::O2>>>
        is_extra_node(_vei_count);
    for (size_t vei=0; vei!=_vei_count; ++vei) {
        for (size_t i=0; i!=VTX_PAIRS_VOL.size(); ++i)
        {
            const std::tuple<size_t,size_t> tup = make_ascending_tuple(
                _vei_to_ni[vei][VTX_PAIRS_VOL[i][0]],
                _vei_to_ni[vei][VTX_PAIRS_VOL[i][1]]
            );
            if (!idxs_extra_nodes.contains(tup)) {
                is_extra_node[vei][i] = true;
                _vei_to_ni[vei][NODES_IN_VOL_ELEM<ElementOrder::O1> + i] =
                    _ni_count;
                idxs_extra_nodes.insert({tup, _ni_count});
                ++_ni_count;
            } else {
                is_extra_node[vei][i] = false;          
                _vei_to_ni[vei][NODES_IN_VOL_ELEM<ElementOrder::O1> + i] =
                    idxs_extra_nodes.at(tup);
            }     
        }
    }

    // TODO: grow() here
    auto temp = std::move(_ni_to_coords);
    _ni_to_coords = fz::SafePtr<std::array<double,3>>(_ni_count);
    std::copy(temp.begin(), temp.end(), _ni_to_coords.begin());
    temp.free();
    
    // second pass: create the extra nodes and assign to volume elements
    for (size_t vei=0; vei!=_vei_count; ++vei) {
        for (size_t i=0; i!=VTX_PAIRS_VOL.size(); ++i)
        {   
            if (!is_extra_node[vei][i]) { continue; }

            const std::tuple<size_t,size_t> tup = make_ascending_tuple(
                _vei_to_ni[vei][VTX_PAIRS_VOL[i][0]],
                _vei_to_ni[vei][VTX_PAIRS_VOL[i][1]]
            );
            const double x = mean(
                _ni_to_coords[std::get<0>(tup)][0],
                _ni_to_coords[std::get<1>(tup)][0]
            );
            const double y = mean(
                _ni_to_coords[std::get<0>(tup)][1],
                _ni_to_coords[std::get<1>(tup)][1]
            );
            const double z = mean(
                _ni_to_coords[std::get<0>(tup)][2],
                _ni_to_coords[std::get<1>(tup)][2]
            );
            const size_t idx_extra_node = idxs_extra_nodes.at(tup);
            _ni_to_coords[idx_extra_node] = {x, y, z}; // add extra node
        }
    }
    is_extra_node.free();

    // third pass: assign nodes to surface elements
    for (size_t sei=0; sei!=_sei_count; ++sei) {
        for (size_t i=0; i!=VTX_PAIRS_SFC.size(); ++i)
        {
            // create a tuple of the indices in ascending order
            const std::tuple<size_t,size_t> tup = make_ascending_tuple(
                _sei_to_ni[sei][VTX_PAIRS_SFC[i][0]],
                _sei_to_ni[sei][VTX_PAIRS_SFC[i][1]]
            );
            _sei_to_ni[sei][NODES_IN_SFC_ELEM<ElementOrder::O1> + i] =
                idxs_extra_nodes.at(tup);
        }
    }
}

// explicit instantiation declarations
template class Simulation<
    Phenomenon::ACOUSTIC,
    NumericalMethod::FEM,
    Domain::FREQUENCY,
    Dimension::D3,
    ElementOrder::O1
>;
template class Simulation<
    Phenomenon::ACOUSTIC,
    NumericalMethod::FEM,
    Domain::FREQUENCY,
    Dimension::D3,
    ElementOrder::O2
>;

} // namespace numav
