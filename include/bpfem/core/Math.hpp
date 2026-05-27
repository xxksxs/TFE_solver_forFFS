#pragma once

#include "bpfem/core/Types.hpp"

#include <array>

namespace fem {

Vec3 operator+(const Vec3& a, const Vec3& b);
Vec3 operator-(const Vec3& a, const Vec3& b);
Vec3 operator*(double s, const Vec3& a);
Vec3 operator/(const Vec3& a, double s);
double dot(const Vec3& a, const Vec3& b);
Vec3 cross(const Vec3& a, const Vec3& b);
double norm(const Vec3& a);
double component(const Vec3& a, int axis);
double tetraVolume(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d);
std::array<Vec3, 4> tetraGradients(const std::array<Vec3, 4>& p);

}  // namespace fem
