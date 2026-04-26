// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/nmvr-format.hpp"

namespace numav::nmvr {

void write_data_chunk(
    std::ofstream& file,
    const std::array<char,CHUNK_ID_SIZE> chunk_id,
    const uint64_t chunk_size_in_bytes,
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

void write_data_chunk_header(
    std::ofstream& file,
    const std::array<char,CHUNK_ID_SIZE> chunk_id,
    const uint64_t chunk_size_in_bytes
) {
    file.write(
        chunk_id.data(), CHUNK_ID_SIZE
    );
    file.write(
        reinterpret_cast<const char*>(&chunk_size_in_bytes), sizeof(uint64_t)
    );
}

void write_data_chunk_body(
    std::ofstream& file,
    const uint64_t size_in_bytes_to_write,
    const void* const data
) {
    file.write(
        reinterpret_cast<const char*>(data), size_in_bytes_to_write
    );
}

} // namespace numav::nmvr
