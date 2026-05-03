// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "common/maths.hpp"
#include "common/utils.hpp"

#include <cassert>

#include "Eigen/Eigen"

namespace numav {

fz::SafePtr<Float> linspace(
    const Float start, const Float finish, const size_t num_points
) {
    assert(num_points != 0UL && num_points != 1UL);
    fz::SafePtr<Float> result(num_points);
    Float step = (finish - start) / static_cast<Float>(num_points - 1UL);
    result.front() = start;
    for (Float* it = result.begin()+1UL; it != result.end()-1UL; ++it) {
        *it = *(it - 1UL) + step;
    }
    result.back() = finish;
    return result;
}

fz::SafePtr<Float> cubspace(
    const Float start, const Float finish, const size_t num_points
) {
    assert(num_points != 0UL && num_points != 1UL);
    fz::SafePtr<Float> x = linspace(0_F, 1_F, num_points);
    const Float& a3 = power<3UL>(start);
    const Float& b3 = power<3UL>(finish);
    fz::SafePtr<Float> result(num_points);
    for (size_t i = 0UL; i != num_points; ++i) {
        result[i] = std::pow(a3 + (b3 - a3) * x[i], 1_F/3_F);
    }
    x.free();
    return result;
}

Float get_triangle_area(const std::array<std::array<Float,3UL>,3UL> coords)
{
    const Float a = std::sqrt(
        power<2UL>(coords[0UL][0UL] - coords[1UL][0UL]) +
        power<2UL>(coords[0UL][1UL] - coords[1UL][1UL]) +
        power<2UL>(coords[0UL][2UL] - coords[1UL][2UL])
    );
    const Float b = std::sqrt(
        power<2UL>(coords[1UL][0UL] - coords[2UL][0UL]) +
        power<2UL>(coords[1UL][1UL] - coords[2UL][1UL]) +
        power<2UL>(coords[1UL][2UL] - coords[2UL][2UL])
    );
    const Float c = std::sqrt(
        power<2UL>(coords[2UL][0UL] - coords[0UL][0UL]) +
        power<2UL>(coords[2UL][1UL] - coords[0UL][1UL]) +
        power<2UL>(coords[2UL][2UL] - coords[0UL][2UL])
    );
    const Float s = (a + b + c) / 2_F;
    const Float area = std::sqrt(s * (s-a) * (s-b) * (s-c));
    return area;
}

Float get_tetrahedron_volume(
    const std::array<std::array<Float,3UL>,4UL> coords
) {
    const Eigen::Matrix<Float,4UL,4UL> mat {
      {coords[0UL][0UL], coords[1UL][0UL], coords[2UL][0UL], coords[3UL][0UL]},
      {coords[0UL][1UL], coords[1UL][1UL], coords[2UL][1UL], coords[3UL][1UL]},
      {coords[0UL][2UL], coords[1UL][2UL], coords[2UL][2UL], coords[3UL][2UL]},
      {             1_F,              1_F,              1_F,              1_F}
    };
    const Float volume = (1_F/6_F) * std::abs(mat.determinant());
    return volume;
}

} // namespace numav
