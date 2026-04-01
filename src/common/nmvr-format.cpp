// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/nmvr-format.hpp"

namespace numav::nmvr {

void write_data_chunk(
    std::ofstream& file,
    const std::array<char,CHUNK_ID_SIZE>& chunk_id,
    const uint64_t& chunk_size_in_bytes,
    const void* const chunk_data
) {
    file.write(
        chunk_id.data(), CHUNK_ID_SIZE
    );
    file.write(
        reinterpret_cast<const char*>(&chunk_size_in_bytes), sizeof(uint64_t)
    );
    file.write(
        reinterpret_cast<const char*>(chunk_data), chunk_size_in_bytes
    );
}

} // namespace numav::nmvr
