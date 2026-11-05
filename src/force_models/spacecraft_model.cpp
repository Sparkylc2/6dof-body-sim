#include "force_models/spacecraft_model.hpp"
#include <algorithm>
#include <cmath>

namespace BSIM {

// --- mass properties ---
MassProperties mass_properties(const SpacecraftParams &p, double m) {
    // incase we
    const double m_fuel = m - p.m_empty;

    MassProperties mp;
    mp.m = p.m_empty + m_fuel;

    // the body origin is the empty cm, so the cm will only move cuz of fuel loss
    mp.r_cm = m_fuel / mp.m * p.r_tank;

    // parallel axis
    auto shift = [](double mass, const Vector3d &d) {
        return mass * (d.squaredNorm() * Matrix3d::Identity() - d * d.transpose());
    };
    // i'm assuming fuel is a solid sphere, which gives J = 2/5 m R^2
    const Matrix3d J_fuel = 0.4 * m_fuel * p.R_tank * p.R_tank * Matrix3d::Identity();

    mp.J = p.J_empty + shift(p.m_empty, -mp.r_cm) + J_fuel + shift(m_fuel, p.r_tank - mp.r_cm);

    // the wheels' spin axis inertia is accounted for in l, so it must be taken out of J
    // ie J = J_spacecraft + J_wheels
    // sum of J_w*a*a^T over the wheels is J_w*A*A^T
    const Matrix3d &A = p.wheels.A;
    mp.J_s = mp.J - p.wheels.J_w * A * A.transpose();

    return mp;
}

// --- gravity ---
Vector3d two_body_gravity(const Vector3d &r) {
    const double rn = r.norm();
    return -MU_EARTH / (rn * rn * rn) * r;
}

Vector3d j2_gravity(const Vector3d &r) {
    const double rn = r.norm();
    const double z2 = r.z() * r.z() / (rn * rn); // (z/r)^2

    // -grad of the J2 term in the potential
    // k is just the fat coefficient at the front
    const double k = -1.5 * J2_EARTH * MU_EARTH * R_EARTH * R_EARTH / std::pow(rn, 5);
    return k *
           Vector3d(r.x() * (1.0 - 5.0 * z2), r.y() * (1.0 - 5.0 * z2), r.z() * (3.0 - 5.0 * z2));
}

double gravity_potential(const Vector3d &r, bool j2) {
    const double rn = r.norm();
    double U = -MU_EARTH / rn;
    if (j2) {
        U += 0.5 * MU_EARTH * J2_EARTH * R_EARTH * R_EARTH *
             (3.0 * r.z() * r.z() / std::pow(rn, 5) - 1.0 / std::pow(rn, 3));
    }
    return U;
}

Vector3d gravity_gradient_torque(const Vector3d &r_b, const Matrix3d &J) {
    // 3*mu / |r|^5 (r cross J*r)
    return 3.0 * MU_EARTH / std::pow(r_b.norm(), 5) * r_b.cross(J * r_b);
}

// --- actuators ---
ThrusterOutput thrust(const SpacecraftParams &p, const std::vector<double> &throttle,
                      const Vector3d &r_cm, double m) {
    ThrusterOutput out;

    // if the tanks empty the thrusters shouldn't work
    if (m <= p.m_empty) {
        return out;
    }

    for (size_t k = 0; k < p.thrusters.size() && k < throttle.size(); ++k) {
        const Thruster &t = p.thrusters[k];
        const double T = std::clamp(throttle[k], 0.0, 1.0) * t.T_max;
        const Vector3d F = T * t.dir;

        out.wrench.F += F;
        // the torque is about wherever the cm is right now, not the body origin
        out.wrench.M += (t.pos - r_cm).cross(F);
        // T = mdot * v_exhaust and v_exhaust = Isp g0
        out.m_dot -= T / (t.Isp * G0);
    }

    return out;
}

Vector3d wheel_motor_torque(const ReactionWheels &w, const Vector3d &tau_cmd, const Vector3d &l) {
    Vector3d tau = tau_cmd.cwiseMax(-w.tau_max).cwiseMin(w.tau_max);

    // if a wheels at top speed you can't push it further so return 0
    for (int k = 0; k < 3; ++k) {
        if ((l(k) >= w.h_max && tau(k) > 0.0) || (l(k) <= -w.h_max && tau(k) < 0.0)) {
            tau(k) = 0.0;
        }
    }
    return tau;
}

Vector3d wheel_speeds(const ReactionWheels &w, const Vector3d &l, const Vector3d &omega) {
    // l_k = J_w*(a_k.w + Omega_k)
    return l / w.J_w - w.A.transpose() * omega;
}

// --- orbits ---
std::pair<Vector3d, Vector3d> circular_orbit(double altitude, double inclination) {
    const double r = R_EARTH + altitude;
    // gravity is the centripetal force, mu/r^2 = v^2 / r
    const double v = std::sqrt(MU_EARTH / r);
    return {Vector3d(r, 0.0, 0.0), v * Vector3d(0.0, std::cos(inclination), std::sin(inclination))};
}

} // namespace BSIM
