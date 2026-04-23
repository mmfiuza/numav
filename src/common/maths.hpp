// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "common/aliases.hpp"

#include "SafePtr.hpp"

#include <tuple>

namespace numav {

fz::SafePtr<double> linspace(
    const double& start, const double& finish, const size_t& num_points
);

fz::SafePtr<double> cubspace(
    const double& start, const double& finish, const size_t& num_points
);

double get_triangle_area(
    const std::array<std::array<double,3UL>,3UL>& coords
);

double get_tetrahedron_volume(
    const std::array<std::array<double,3UL>,4UL>& coords
);

template<typename... T>
auto mean(const T&... args) {
    return (args + ...) / (sizeof...(args));
}

template<size_t N> constexpr size_t FACTORIAL = [] {
    size_t result = 1UL;
    for(size_t i = 2UL; i <= N; ++i) {
        result *= i;
    }
    return result;
}();

template<size_t N, size_t K> constexpr size_t COMB_REP_SIZE = [] {
    // combination with repetition
    return FACTORIAL<N+K-1UL> / (FACTORIAL<K> * FACTORIAL<N-1UL>);
}();

template<size_t N>
constexpr std::array<std::array<size_t,2UL>, COMB_REP_SIZE<N,2UL>> 
COMBINATION_REP = [] { // todo: generalize for K!=2 (maybe)
    constexpr size_t K = 2UL;
    // combination with repetition in upper column major order
    std::array<std::array<size_t,K>, COMB_REP_SIZE<N,K>> result;
    auto it_result = result.begin();
    for (size_t j = 0UL; j != N; ++j) {
        for (size_t i = 0UL; i != j+1UL; ++i) {
            *it_result = {i,j};
            ++it_result;
        }
    }
    return result;
}();

template<typename XContainer, typename YContainer, typename XType>
auto interpolate(
    const XContainer& x_array,
    const YContainer& y_array,
    const XType& x_value)
{
    // Obtain sizes and iterators in a container‑agnostic way
    auto x_size = std::size(x_array);
    auto y_size = std::size(y_array);
    auto x_begin = std::begin(x_array);
    auto x_end   = std::end(x_array);
    auto y_begin = std::begin(y_array);

    // Input validation
    if (x_size == 0UL || y_size == 0UL)
        throw std::invalid_argument("Vectors must not be empty.");
    if (x_size != y_size)
        throw std::invalid_argument(
            "x_array and y_array must have the same size."
        );

    // Check strict monotonicity of x_array
    if (x_size >= 2UL) {
        auto prev = x_begin;
        auto curr = std::next(prev);
        while (curr != x_end) {
            if (*prev >= *curr)
                throw std::invalid_argument(
                    "x_array must be strictly increasing."
                );
            ++prev;
            ++curr;
        }
    }

    // Single point case
    if (x_size == 1UL) {
        if (x_value == *x_begin)
            return *y_begin;
        else
            throw std::out_of_range("Extrapolation not allowed.");
    }

    // Reject extrapolation
    if (x_value < *x_begin || x_value > *std::prev(x_end))
        throw std::out_of_range("Extrapolation not allowed.");

    // Locate the interval using binary search
    auto it = std::lower_bound(x_begin, x_end, x_value);

    // Exact match at a knot (including boundaries)
    if (it != x_end && x_value == *it) {
        auto index = std::distance(x_begin, it);
        auto y_it = y_begin;
        std::advance(y_it, index);
        return *y_it;
    }

    // Otherwise, interval is between the previous element and the found one
    // Note: it cannot be x_begin because we already handled x_value < *x_begin
    auto prev_it = std::prev(it); // guaranteed to exist because it != x_begin

    // Get iterators to the corresponding y values
    auto prev_index = std::distance(x_begin, prev_it);
    auto curr_index = std::distance(x_begin, it);
    auto y_prev = y_begin;
    auto y_curr = y_begin;
    std::advance(y_prev, prev_index);
    std::advance(y_curr, curr_index);

    // Linear interpolation formula
    return *y_prev + (*y_curr - *y_prev) *
           (x_value - *prev_it) / (*it - *prev_it);
}

} // namespace numav
