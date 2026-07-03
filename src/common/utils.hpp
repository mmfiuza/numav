// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/numav.hpp"

#include <string_view>
#include <charconv>
#include <tuple>
#include <fstream>
#include <memory>

namespace numav {

// TODO: check if pass by ref is needed
void trim_right_whitespace(std::string_view& sv); 

template<typename T>
T parse(std::string_view str) {
    trim_right_whitespace(str);
    const std::string tmp(str); // strto* requires null-termination
    char* end;
    errno = 0;
    T value;
    if constexpr (std::is_same_v<T, float>) {
        value = std::strtof(tmp.c_str(), &end);
    }
    else if constexpr (std::is_same_v<T, double>) {
        value = std::strtod(tmp.c_str(), &end);
    }
    else if constexpr (std::is_same_v<T, long double>) {
        value = std::strtold(tmp.c_str(), &end);
    }
    else if constexpr (std::is_integral_v<T> && std::is_signed_v<T>) {
        value = static_cast<T>(std::strtoll(tmp.c_str(), &end, 10));
    }
    else if constexpr (std::is_integral_v<T> && std::is_unsigned_v<T>) {
        value = static_cast<T>(std::strtoull(tmp.c_str(), &end, 10));
    }
    else {
        static_assert(!sizeof(T), "unsupported type for parse()");
    }
    if (end != tmp.c_str() + tmp.size() || errno == ERANGE) {
        throw std::invalid_argument("invalid number format");
    }
    return value;
}

FuncFloatToCmplx const2func(
    const Cmplx constant
);

FuncFloatToCmplx table2func(
    const char* const impedance_text_file
);

template<typename T>
std::tuple<T,T> make_ascending_tuple(const T a, const T b) {
    return a<b ? std::make_tuple(a,b) : std::make_tuple(b,a);
}

template<typename T, std::uint64_t... Sizes>
constexpr auto concat_constexpr_arrays(const std::array<T, Sizes>... arrays)
{
    constexpr std::uint64_t total_size = (Sizes + ...);
    std::array<T, total_size> result{};
    std::uint64_t index = 0UL;
    (
        (
            [&]() {
                for (std::uint64_t i = 0UL; i < arrays.size(); ++i) {
                    result[index + i] = arrays[i];
                }
            }(),
            index += arrays.size()
        ), ...
    );
    return result;
}

template<typename T1, typename T2>
std::unique_ptr<T1[]> static_cast_contiguous_data(
    const T2* const data,
    const uint64_t element_count
) {
    std::unique_ptr<T1[]> ptr = std::make_unique<T1[]>(element_count);
    for (uint64_t i = 0UL; i!= element_count; ++i){
        ptr[i] = static_cast<T1>(data[i]);
    }
    return ptr;
}

template< typename T>
std::unique_ptr<uint64_t[]> to_one_based_index(
    const T* const data,
    const uint64_t element_count
) {
    std::unique_ptr<uint64_t[]> ptr = std::make_unique<uint64_t[]>(element_count);
    for (uint64_t i = 0UL; i!= element_count; ++i){
        ptr[i] = data[i] + 1UL;
    }
    return ptr;
}

} // namespace numav
