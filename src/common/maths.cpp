// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/maths.hpp"

#include <cassert>

#include "Eigen/Eigen"

namespace numav {

fz::SafePtr<double> linspace(
    const double& start, const double& finish, const size_t& num_points
) {
    assert(num_points != 0UL && num_points != 1UL);
    fz::SafePtr<double> result(num_points);
    double step = (finish - start) / static_cast<double>(num_points - 1UL);
    result.front() = start;
    for (double* it=result.begin()+1UL; it!=result.end()-1UL; ++it) {
        *it = *(it - 1UL) + step;
    }
    result.back() = finish;
    return result;
}

fz::SafePtr<double> cubspace(
    const double& start, const double& finish, const size_t& num_points
) {
    assert(num_points!=0 && num_points!=1);
    fz::SafePtr<double> x = linspace(0.0, 1.0, num_points);
    const double& a3 = std::pow(start, 3UL);
    const double& b3 = std::pow(finish, 3UL);
    fz::SafePtr<double> result(num_points);
    for (size_t i = 0UL; i != num_points; ++i) {
        result[i] = std::pow(a3 + (b3 - a3) * x[i], 1.0/3.0);
    }
    x.free();
    return result;
}

double get_triangle_area(const std::array<std::array<double,3UL>,3UL>& coords)
{
    const double a = std::sqrt(
        std::pow(coords[0UL][0UL] - coords[1UL][0UL], 2UL) +
        std::pow(coords[0UL][1UL] - coords[1UL][1UL], 2UL) +
        std::pow(coords[0UL][2UL] - coords[1UL][2UL], 2UL)
    );
    const double b = std::sqrt(
        std::pow(coords[1UL][0UL] - coords[2UL][0UL], 2UL) +
        std::pow(coords[1UL][1UL] - coords[2UL][1UL], 2UL) +
        std::pow(coords[1UL][2UL] - coords[2UL][2UL], 2UL)
    );
    const double c = std::sqrt(
        std::pow(coords[2UL][0UL] - coords[0UL][0UL], 2UL) +
        std::pow(coords[2UL][1UL] - coords[0UL][1UL], 2UL) +
        std::pow(coords[2UL][2UL] - coords[0UL][2UL], 2UL)
    );
    const double s = (a + b + c) / 2.0;
    const double area = std::sqrt(s * (s-a) * (s-b) * (s-c));
    return area;
}

double get_tetrahedron_volume(
    const std::array<std::array<double,3UL>,4UL>& coords
) {
    const Eigen::Matrix<double,4UL,4UL> mat {
      {coords[0UL][0UL], coords[1UL][0UL], coords[2UL][0UL], coords[3UL][0UL]},
      {coords[0UL][1UL], coords[1UL][1UL], coords[2UL][1UL], coords[3UL][1UL]},
      {coords[0UL][2UL], coords[1UL][2UL], coords[2UL][2UL], coords[3UL][2UL]},
      {             1.0,              1.0,              1.0,              1.0}
    };
    const double volume = (1.0/6.0) * std::abs(mat.determinant());
    return volume;
}

} // namespace numav
