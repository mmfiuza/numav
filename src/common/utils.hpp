// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include <string_view>
#include <charconv>
#include <stdexcept>
#include <tuple>
#include <fstream>

void trim_right_whitespace(std::string_view& sv);

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

template<typename T>
std::tuple<T,T> make_ascending_tuple(const T& a, const T& b) {
    return a<b ? std::make_tuple(a,b) : std::make_tuple(b,a);
}

template<typename T>
std::pair<T,T> make_ascending_pair(const T& a, const T& b) {
    return a<b ? std::make_pair(a,b) : std::make_pair(b,a);
}

template<typename T>
bool compare_pair(const std::pair<T,T>& a, const std::pair<T,T>& b) {
    #if NUMAV_GLOBAL_MATRIX_STORAGE_ORDER == NUMAV_UPPER_COL_MAJOR
        if (a.second != b.second) {
            return a.second < b.second;
        }
        else {
            return a.first < b.first;
        }
    #elif NUMAV_GLOBAL_MATRIX_STORAGE_ORDER == NUMAV_UPPER_ROW_MAJOR
        if (a.first != b.first) {
            return a.first < b.first;
        }
        else {
            return a.second < b.second;
        }
    #else
        static_assert(false, "Invalid GLOBAL_MATRIX_STORAGE_ORDER.");
    #endif
}

template<typename T>
void write_matrix(const T& matrix, const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (file.is_open()) {
        const uint64_t rows = matrix.rows();
        const uint64_t cols = matrix.cols();
        file.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
        file.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
        file.write(
            reinterpret_cast<const char*>(matrix.data()),
            rows * cols * sizeof(typename T::Scalar)
        );
    }
}
