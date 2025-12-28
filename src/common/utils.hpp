// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "common/debug-macros.hpp"
#include "common/aliases.hpp"

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

template<typename T, std::size_t... Sizes>
constexpr auto concat_constexpr_arrays(const std::array<T, Sizes>&... arrays)
{
    constexpr std::size_t total_size = (Sizes + ...);
    std::array<T, total_size> result{};
    
    std::size_t index = 0;
    (
        (
            std::copy(arrays.begin(), arrays.end(), result.begin() + index), 
            index += arrays.size()
        ), ...
    );
    
    return result;
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
