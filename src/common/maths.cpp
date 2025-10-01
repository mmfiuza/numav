// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/maths.hpp"

#include <cassert>

fz::SafePtr<double> linspace(
    const double& start, const double& finish, const size_t& num_points
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

double get_triangle_area(const std::array<std::array<double,3>,3>& coords)
{
    const double a = std::sqrt(
        std::pow(coords[0][0] - coords[1][0], 2) +
        std::pow(coords[0][1] - coords[1][1], 2) +
        std::pow(coords[0][2] - coords[1][2], 2)
    );
    const double b = std::sqrt(
        std::pow(coords[1][0] - coords[2][0], 2) +
        std::pow(coords[1][1] - coords[2][1], 2) +
        std::pow(coords[1][2] - coords[2][2], 2)
    );
    const double c = std::sqrt(
        std::pow(coords[2][0] - coords[0][0], 2) +
        std::pow(coords[2][1] - coords[0][1], 2) +
        std::pow(coords[2][2] - coords[0][2], 2)
    );
    const double s = (a + b + c) / 2.0;
    const double area = std::sqrt(s * (s-a) * (s-b) * (s-c));
    return area;
}
