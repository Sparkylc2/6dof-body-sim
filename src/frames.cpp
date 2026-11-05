#include "frames.hpp"
#include <algorithm>
#include <cmath>

namespace BSIM {

Matrix3d rotation_euler(double phi, double theta, double psi) {
    const double c_ph = std::cos(phi), s_ph = std::sin(phi);
    const double c_th = std::cos(theta), s_th = std::sin(theta);
    const double c_ps = std::cos(psi), s_ps = std::sin(psi);

    // this comes from multiplying the three rotation matrices
    // for the three different rotations (eqn 2.5)
    const Matrix3d R{
        {c_th * c_ps, c_th * s_ps, -s_th},
        {s_ph * s_th * c_ps - c_ph * s_ps, s_ph * s_th * s_ps + c_ph * c_ps, s_ph * c_th},
        {c_ph * s_th * c_ps + s_ph * s_ps, c_ph * s_th * s_ps - s_ph * c_ps, c_ph * c_th},
    };

    return R;
}

Matrix3d body_rates_to_euler_rates(double phi, double theta) {
    const double s_phi = std::sin(phi), c_phi = std::cos(phi);
    const double t_theta = std::tan(theta), c_theta = std::cos(theta);

    // exactly whats in eqn 3.16, from differentiating and the rotation procedure
    const Matrix3d R{
        {1.0, s_phi * t_theta, c_phi * t_theta},
        {0.0, c_phi, -s_phi},
        {0.0, s_phi / c_theta, c_phi / c_theta},
    };
    return R;
}

Matrix3d rotation(const Vector4d &e) {
    const double e0 = e[0], e1 = e[1], e2 = e[2], e3 = e[3];

    // verbatim from B.1, this is R_b^v (body to vehicle/inertial)
    const Matrix3d R{
        {e1 * e1 + e0 * e0 - e2 * e2 - e3 * e3, 2 * (e1 * e2 - e3 * e0), 2 * (e1 * e3 + e2 * e0)},
        {2 * (e1 * e2 + e3 * e0), e2 * e2 + e0 * e0 - e1 * e1 - e3 * e3, 2 * (e2 * e3 - e1 * e0)},
        {2 * (e1 * e3 - e2 * e0), 2 * (e2 * e3 + e1 * e0), e3 * e3 + e0 * e0 - e1 * e1 - e2 * e2},
    };

    return R;
}

Vector4d quat_rate(const Vector4d &e, const Vector3d &w) {
    const double p = w.x(), q = w.y(), r = w.z();

    // verbatim from B.3
    const Eigen::Matrix4d Omega{
        {0, -p, -q, -r},
        {p, 0, r, -q},
        {q, -r, 0, p},
        {r, q, -p, 0},
    };

    // e_dot
    return 0.5 * Omega * e;
}

Vector4d euler_to_quat(double phi, double theta, double psi) {
    const double c_phi = std::cos(phi / 2), s_phi = std::sin(phi / 2);
    const double c_theta = std::cos(theta / 2), s_theta = std::sin(theta / 2);
    const double c_psi = std::cos(psi / 2), s_psi = std::sin(psi / 2);

    // composing the three half angle rotations (psi about z, then theta about y, then phi about x)
    return Vector4d{
        c_psi * c_theta * c_phi + s_psi * s_theta * s_phi,
        c_psi * c_theta * s_phi - s_psi * s_theta * c_phi,
        c_psi * s_theta * c_phi + s_psi * c_theta * s_phi,
        s_psi * c_theta * c_phi - c_psi * s_theta * s_phi,
    };
}

Vector3d quat_to_euler(const Vector4d &e) {
    const double e0 = e[0], e1 = e[1], e2 = e[2], e3 = e[3];

    // matching entries of B.1 against eqn 2.5
    const double phi = std::atan2(2 * (e0 * e1 + e2 * e3), e0 * e0 + e3 * e3 - e1 * e1 - e2 * e2);
    // clamp so rounding near +-90 deg pitch isn't able to push asin outside [-1, 1]
    const double theta = std::asin(std::clamp(2 * (e0 * e2 - e1 * e3), -1.0, 1.0));
    const double psi = std::atan2(2 * (e0 * e3 + e1 * e2), e0 * e0 + e1 * e1 - e2 * e2 - e3 * e3);

    return Vector3d{phi, theta, psi};
}

} // namespace BSIM
