// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/maths.hpp"

#include <cassert>

#include "Eigen/Eigen"

namespace numav {

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

fz::SafePtr<double> cubespace(
    const double& start, const double& finish, const size_t& num_points
) {
    assert(num_points!=0 && num_points!=1);
    fz::SafePtr<double> x = linspace(0.0, 1.0, num_points);
    const double& a3 = std::pow(start, 3);
    const double& b3 = std::pow(finish, 3);
    fz::SafePtr<double> result(num_points);
    for (size_t i=0; i!=num_points; ++i) {
        result[i] = std::pow(a3 + (b3 - a3) * x[i], 1.0/3.0);
    }
    x.free();
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

double get_tetrahedron_volume(const std::array<std::array<double,3>,4>& coords)
{
    const Eigen::Matrix<double,4,4> mat {
        {coords[0][0], coords[1][0], coords[2][0], coords[3][0]},
        {coords[0][1], coords[1][1], coords[2][1], coords[3][1]},
        {coords[0][2], coords[1][2], coords[2][2], coords[3][2]},
        {           1,            1,            1,            1}
    };
    const double volume = (1.0/6.0) * std::abs(mat.determinant());
    return volume;
}

} // namespace numav
