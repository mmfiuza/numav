// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav/numav.hpp"

#include <tuple>
#include <fstream>
#include <limits>

#include "common/log.hpp"
#include "common/hash-functions.hpp"
#include "common/maths.hpp"
#include "common/utils.hpp"

namespace numav {

template <ElementOrder O>
_idx_t SimulationAcFemFreqD3<O>::_get_closest_point(
    const std::array<double,3>& point_coords
) {
    double minimum_distance_squared = std::numeric_limits<double>::max();
    _idx_t ni_closest = std::numeric_limits<size_t>::max();
    for (_idx_t ni=0; ni!=_ni_count(); ++ni) {
        double distance_squared = 
            std::pow(_node_coords[ni][0] - point_coords[0], 2) +
            std::pow(_node_coords[ni][1] - point_coords[1], 2) +
            std::pow(_node_coords[ni][2] - point_coords[2], 2)
        ;
        if (distance_squared < minimum_distance_squared) {
            minimum_distance_squared = distance_squared;
            ni_closest = ni;
        }
    }
    std::cout << "ni_closest: " << ni_closest << "\n";
    return ni_closest;
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
        fz::SafePtr<std::array<size_t,NODES_IN_SFC_ELEM<O>>>(sfc_elem_count);
    _vol_elem_node_idx =
        fz::SafePtr<std::array<size_t,NODES_IN_VOL_ELEM<O>>>(vol_elem_count);
    _sei_to_espg = fz::SafePtr<size_t>(sfc_elem_count);
    _vei_to_evpg = fz::SafePtr<size_t>(vol_elem_count);
    auto it_node_coords = _node_coords.begin();
    auto it_sfc_elem_node_idx = _sfc_elem_node_idx.begin();
    auto it_vol_elem_node_idx = _vol_elem_node_idx.begin();
    auto it_elem_idx_to_espg = _sei_to_espg.begin();
    auto it_elem_idx_to_evpg = _vei_to_evpg.begin();
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
            const _epg_t espg = parse<size_t>(line.substr(16,8));
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
            const _epg_t evpg = parse<size_t>(line.substr(16,8));
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
void SimulationAcFemFreqD3<ElementOrder::O1>::_generate_extra_nodes() {
    // nothing needs to be done in this case (in order 1)
}

template<>
void SimulationAcFemFreqD3<ElementOrder::O2>::_generate_extra_nodes()
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
        is_extra_node(_vei_count());
    size_t count = _ni_count();
    for (size_t e=0; e!=_vei_count(); ++e) {
        for (size_t i=0; i!=VTX_PAIRS_VOL.size(); ++i)
        {
            const std::tuple<size_t,size_t> tup = make_ascending_tuple(
                _vol_elem_node_idx[e][VTX_PAIRS_VOL[i][0]],
                _vol_elem_node_idx[e][VTX_PAIRS_VOL[i][1]]
            );
            if (!idxs_extra_nodes.contains(tup)) {
                is_extra_node[e][i] = true;
                _vol_elem_node_idx[e][NODES_IN_VOL_ELEM<ElementOrder::O1> + i] =
                    count;
                idxs_extra_nodes.insert({tup, count});
                ++count;
            } else {
                is_extra_node[e][i] = false;          
                _vol_elem_node_idx[e][NODES_IN_VOL_ELEM<ElementOrder::O1> + i] =
                    idxs_extra_nodes.at(tup);
            }     
        }
    }

    // TODO: grow() here
    auto temp = std::move(_node_coords);
    _node_coords = fz::SafePtr<std::array<double,3>>(count);
    std::copy(temp.begin(), temp.end(), _node_coords.begin());
    temp.free();
    
    // second pass: create the extra nodes and assign to volume elements
    for (size_t e=0; e!=_vei_count(); ++e) {
        for (size_t i=0; i!=VTX_PAIRS_VOL.size(); ++i)
        {   
            if (!is_extra_node[e][i]) { continue; }

            const std::tuple<size_t,size_t> tup = make_ascending_tuple(
                _vol_elem_node_idx[e][VTX_PAIRS_VOL[i][0]],
                _vol_elem_node_idx[e][VTX_PAIRS_VOL[i][1]]
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

    // third pass: assign nodes to surface elements
    for (size_t e=0; e!=_sei_count(); ++e) {
        for (size_t i=0; i!=VTX_PAIRS_SFC.size(); ++i)
        {
            // create a tuple of the indices in ascending order
            const std::tuple<size_t,size_t> tup = make_ascending_tuple(
                _sfc_elem_node_idx[e][VTX_PAIRS_SFC[i][0]],
                _sfc_elem_node_idx[e][VTX_PAIRS_SFC[i][1]]
            );
            _sfc_elem_node_idx[e][NODES_IN_SFC_ELEM<ElementOrder::O1> + i] =
                idxs_extra_nodes.at(tup);
        }
    }
}

} // namespace numav