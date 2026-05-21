// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/aliases.hpp"

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
    T value;
    auto result = std::from_chars(str.data(), str.data() + str.size(), value);
    if (result.ec != std::errc{} || result.ptr != str.data()+str.size()) {
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

inline constexpr Float operator"" _F(unsigned long long v) noexcept {
    return static_cast<Float>(v);
}

inline constexpr Float operator"" _F(long double v) noexcept {
    return static_cast<Float>(v);
}

template<typename T>
std::tuple<T,T> make_ascending_tuple(const T a, const T b) {
    return a<b ? std::make_tuple(a,b) : std::make_tuple(b,a);
}

template<typename T, std::size_t... Sizes>
constexpr auto concat_constexpr_arrays(const std::array<T, Sizes>... arrays)
{
    constexpr std::size_t total_size = (Sizes + ...);
    std::array<T, total_size> result{};
    
    std::size_t index = 0UL;
    (
        (
            std::copy(arrays.begin(), arrays.end(), result.begin() + index), 
            index += arrays.size()
        ), ...
    );
    
    return result;
}

template<typename T1, typename T2>
std::unique_ptr<T1[]> static_cast_contiguous_data(
    const T2* const data,
    const size_t element_count
) {
    std::unique_ptr<T1[]> ptr = std::make_unique<T1[]>(element_count);
    for (size_t i = 0UL; i!= element_count; ++i){
        ptr[i] = static_cast<T1>(data[i]);
    }
    return ptr;
}

} // namespace numav
