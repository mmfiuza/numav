// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

fz::SafePtr<double> linspace(
    const double& start,
    const double& finish,
    const size_t& num_points
) {
    assert(num_points!=0 && num_points!=1);
    fz::SafePtr<double> result(num_points);
    double step = (finish - start) / static_cast<double>(num_points - 1);
    result.front() = start;
    for (double* it=result.begin()+1; it!=result.end()-1; ++it) {
        *it = *(it-1) + step;
    }
    result.back() = finish;
    return result;
}

template<typename... T>
auto mean(const T&... args) {
    return (args + ...) / (sizeof...(args));
}

template<typename T>
std::tuple<T,T> make_ascending_tuple(const T& a, const T& b) {
    return a<b ? std::make_tuple(a,b) : std::make_tuple(b,a);
}

template<size_t N> constexpr size_t FACTORIAL = [] {
    size_t result = 1;
    for(size_t i = 2; i<=N; ++i) {
        result *= i;
    }
    return result;
}();

template<size_t N, size_t K> constexpr size_t COMB_REP_SIZE = [] {
    // combination with repetition
    return FACTORIAL<N+K-1> / (FACTORIAL<K> * FACTORIAL<N-1>);
}();

template<size_t N>
constexpr std::array<std::array<size_t,2>, COMB_REP_SIZE<N,2>> 
COMBINATION_REP = [] { // todo: generalize for K!=2 (maybe)
    constexpr size_t K = 2;
    // combination with repetition
    std::array<std::array<size_t,K>, COMB_REP_SIZE<N,K>> result;
    auto it_result = result.begin();
    for (size_t j=0; j!=N; ++j) {
        for (size_t i=0; i!=j+1; ++i) {
            *it_result = {i,j};
            ++it_result;
        }
    }
    return result;
}();
