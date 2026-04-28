// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include <fstream>
#include <array>

#include "numav/numav.hpp"

namespace numav::nmvr {

inline constexpr uint64_t SIGNATURE_SIZE = 8UL;
inline constexpr std::array<char, SIGNATURE_SIZE> SIGNATURE = {
    'n','u','m','a','v','r','s','t'
};

inline constexpr uint64_t END_OF_FILE_SIZE = 8UL;
inline constexpr std::array<char, END_OF_FILE_SIZE> END_OF_FILE = {
    'n','u','m','a','v','e','o','f'
};

inline constexpr uint64_t CHUNK_ID_SIZE = 8UL;
inline constexpr std::array<char, CHUNK_ID_SIZE>
    SIMULATION_TYPE_CHUNK_ID = {
        's','i','m','_','t','y','p','e'
    };
inline constexpr std::array<char, CHUNK_ID_SIZE>
    SIMULATED_FREQUENCY_STEPS_CHUNK_ID = {
        'f','r','e','q','_','s','t','p'
    };
inline constexpr std::array<char, CHUNK_ID_SIZE>
    NODE_INDEX_TO_COORDINATES_CHUNK_ID = {
        'n','i','2','c','o','o','r','d'
    };
inline constexpr std::array<char, CHUNK_ID_SIZE>
    SURFACE_ELEMENT_INDEX_TO_NODE_INDEX_CHUNK_ID = {
        's','e','i','_','2','_','n','i'
    };
inline constexpr std::array<char, CHUNK_ID_SIZE>
    VOLUME_ELEMENT_INDEX_TO_NODE_INDEX_CHUNK_ID = {
        'v','e','i','_','2','_','n','i'
    };
inline constexpr std::array<char, CHUNK_ID_SIZE>
    COMPLEX_AMPLITUDE_OF_SOUND_PRESSURE_SOLUTION_CHUNK_ID = {
        'c','p','x','_','p','r','e','s'
    };

inline constexpr uint64_t SIM_TYPE_PHENOMENON_SIZE = 8UL;
inline constexpr std::array<char, SIM_TYPE_PHENOMENON_SIZE>
    SIM_TYPE_PHENOMENON = {
        'a','c','o','u','s','t','i','c'
    };

inline constexpr uint64_t SIM_TYPE_NUMERICAL_METHOD_SIZE = 8UL;
inline constexpr std::array<char, SIM_TYPE_NUMERICAL_METHOD_SIZE>
    SIM_TYPE_NUMERICAL_METHOD = {
        'f','e','m','_','_','_','_','_'
    };

inline constexpr uint64_t SIM_TYPE_DOMAIN_SIZE = 8UL;
template<Domain D>
inline constexpr std::array<char, SIM_TYPE_DOMAIN_SIZE> SIM_TYPE_DOMAIN = [] {
    if constexpr (D == Domain::FREQUENCY) {
        return std::array<char, SIM_TYPE_DOMAIN_SIZE>{
            'f','r','e','q','_','d','o','m'
        };
    }
    if constexpr (D == Domain::TIME) {
        return std::array<char, SIM_TYPE_DOMAIN_SIZE>{
            't','i','m','e','_','d','o','m'
        };
    }
}();

inline constexpr uint64_t SIM_TYPE_DIMENSION_SIZE = 8UL;
template<Dimension D>
inline constexpr std::array<char, SIM_TYPE_DIMENSION_SIZE> SIM_TYPE_DIMENSION =
[] {
    if constexpr (D == Dimension::D1) {
        return std::array<char, SIM_TYPE_DIMENSION_SIZE>{
            '1','_','d','i','m','e','n','s'
        };
    }
    if constexpr (D == Dimension::D2) {
        return std::array<char, SIM_TYPE_DIMENSION_SIZE>{
            '2','_','d','i','m','e','n','s'
        };
    }
    if constexpr (D == Dimension::D3) {
        return std::array<char, SIM_TYPE_DIMENSION_SIZE>{
            '3','_','d','i','m','e','n','s'
        };
    }
}();

inline constexpr uint64_t SIM_TYPE_ELEMENT_ORDER_SIZE = 8UL;
template<ElementOrder O>
inline constexpr std::array<char, SIM_TYPE_ELEMENT_ORDER_SIZE>
SIM_TYPE_FEM_ORDER =
[] {
    if constexpr (O == ElementOrder::O1) {
        return std::array<char, SIM_TYPE_ELEMENT_ORDER_SIZE>{
            '1','s','t','_','o','r','d','_'
        };
    }
    if constexpr (O == ElementOrder::O2) {
        return std::array<char, SIM_TYPE_ELEMENT_ORDER_SIZE>{
            '2','n','d','_','o','r','d','_'
        };
    }
}();

void write_data_chunk(
    std::ofstream&,
    const std::array<char,CHUNK_ID_SIZE>,
    const size_t,
    const void* const
);
void write_data_chunk_header(
    std::ofstream&,
    const std::array<char,CHUNK_ID_SIZE>,
    const size_t
);
void write_data_chunk_body(
    std::ofstream&,
    const size_t,
    const void* const
);

} // namespace numav::nmvr
