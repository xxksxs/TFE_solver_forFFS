#include "bpfem/core/Math.hpp"

#include <cmath>

namespace fem {

Vec3 operator+(const Vec3& a, const Vec3& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 operator-(const Vec3& a, const Vec3& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 operator*(double s, const Vec3& a) {
    return {s * a.x, s * a.y, s * a.z};
}

Vec3 operator/(const Vec3& a, double s) {
    return {a.x / s, a.y / s, a.z / s};
}

double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

double norm(const Vec3& a) {
    return std::sqrt(dot(a, a));
}

double component(const Vec3& a, int axis) {
    if (axis == 0) return a.x;
    if (axis == 1) return a.y;
    return a.z;
}

double tetraVolume(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
    const double ax = b.x - a.x;
    const double ay = b.y - a.y;
    const double az = b.z - a.z;
    const double bx = c.x - a.x;
    const double by = c.y - a.y;
    const double bz = c.z - a.z;
    const double cx = d.x - a.x;
    const double cy = d.y - a.y;
    const double cz = d.z - a.z;
    const double det = ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx);
    return std::abs(det) / 6.0;
}

std::array<Vec3, 4> tetraGradients(const std::array<Vec3, 4>& p) {
    std::array<Vec3, 4> grad{};
    const double x1 = p[0].x, y1 = p[0].y, z1 = p[0].z;
    const double x2 = p[1].x, y2 = p[1].y, z2 = p[1].z;
    const double x3 = p[2].x, y3 = p[2].y, z3 = p[2].z;
    const double x4 = p[3].x, y4 = p[3].y, z4 = p[3].z;
    const double sixV = (x2 - x1) * ((y3 - y1) * (z4 - z1) - (z3 - z1) * (y4 - y1))
                      - (y2 - y1) * ((x3 - x1) * (z4 - z1) - (z3 - z1) * (x4 - x1))
                      + (z2 - z1) * ((x3 - x1) * (y4 - y1) - (y3 - y1) * (x4 - x1));
    if (std::abs(sixV) < 1.0e-30) {
        return grad;
    }

    const std::array<std::array<double, 3>, 4> coeff = {{
        {{y3 * z4 - y4 * z3 + y4 * z2 - y2 * z4 + y2 * z3 - y3 * z2,
          x4 * z3 - x3 * z4 + x2 * z4 - x4 * z2 + x3 * z2 - x2 * z3,
          x3 * y4 - x4 * y3 + x4 * y2 - x2 * y4 + x2 * y3 - x3 * y2}},
        {{y4 * z3 - y3 * z4 + y1 * z4 - y4 * z1 + y3 * z1 - y1 * z3,
          x3 * z4 - x4 * z3 + x4 * z1 - x1 * z4 + x1 * z3 - x3 * z1,
          x4 * y3 - x3 * y4 + x1 * y4 - x4 * y1 + x3 * y1 - x1 * y3}},
        {{y2 * z4 - y4 * z2 + y4 * z1 - y1 * z4 + y1 * z2 - y2 * z1,
          x4 * z2 - x2 * z4 + x1 * z4 - x4 * z1 + x2 * z1 - x1 * z2,
          x2 * y4 - x4 * y2 + x4 * y1 - x1 * y4 + x1 * y2 - x2 * y1}},
        {{y3 * z2 - y2 * z3 + y1 * z3 - y3 * z1 + y2 * z1 - y1 * z2,
          x2 * z3 - x3 * z2 + x3 * z1 - x1 * z3 + x1 * z2 - x2 * z1,
          x3 * y2 - x2 * y3 + x1 * y3 - x3 * y1 + x2 * y1 - x1 * y2}}
    }};

    for (int i = 0; i < 4; ++i) {
        grad[i] = {-coeff[static_cast<std::size_t>(i)][0] / sixV,
                   -coeff[static_cast<std::size_t>(i)][1] / sixV,
                   -coeff[static_cast<std::size_t>(i)][2] / sixV};
    }
    return grad;
}

}  // namespace fem
