#include "autopilots/attitude_control.hpp"
#include "frames.hpp"
#include <cmath>

namespace BSIM {

// --- attitude error ---
Vector4d attitude_error_quat(const Vector4d &q_d, const Vector4d &q) {
    return quat_mul(quat_conj(q_d.normalized()), q);
}

Vector4d pointing_error_quat(const Vector4d &q, const Vector3d &b, const Vector3d &d_i) {
    const Vector3d b_n = b.normalized();
    const Vector3d d_b = rotation(q).transpose() * d_i.normalized(); // target in body axes
    const double c = b_n.dot(d_b);                                   // cos phi

    // if pointing directly away, b x d is 0 so there's no axis,
    // but every perpendicular one is equally short
    if (1.0 + c < 1e-6) {
        Vector3d m = b_n.cross(Vector3d::UnitX());
        if (m.norm() < 0.1) {
            m = b_n.cross(Vector3d::UnitY());
        }
        Vector4d q_e;
        q_e << 0.0, m.normalized();
        return q_e;
    }

    // q_e = q_s* = [1 + c, -(b cross d)] / sqrt(2 (1 + c)), and -(b cross d) = d cross b
    const double k = std::sqrt(2.0 * (1.0 + c));
    Vector4d q_e;
    q_e << (1.0 + c) / k, d_b.cross(b_n) / k;
    return q_e;
}

// --- gains ---
AttitudeGains attitude_gains(const Matrix3d &J_s, const AttitudeDesignParams &d) {
    AttitudeGains g;

    // inner loop
    // J w_dot = k_w (w_c - w) is first order with bandwidth k_w / J
    g.k_w = d.w_rate * J_s.diagonal();
    g.k_wi = d.ki_ratio * d.w_rate * g.k_w;

    // with a P inner loop the two together are a PD with k_p = k_w k_theta and k_d = k_w,
    // so w_n^2 = w_rate k_theta and 2 zeta w_n = w_rate, which gives w_rate / k_theta = 4 zeta^2
    g.k_theta = d.w_rate / (4.0 * d.zeta * d.zeta);
    g.wn = std::sqrt(d.w_rate * g.k_theta);
    return g;
}

// --- attitude controller ---
void AttitudeController::design(const SpacecraftParams &p, double m,
                                const AttitudeDesignParams &d) {
    gains = attitude_gains(mass_properties(p, m).J_s, d);
    A_inv = p.wheels.A.inverse();
    tau_max = p.wheels.tau_max;

    for (int i = 0; i < 3; ++i) {
        PIDParams &o = outer[i].params;
        o.kp = gains.k_theta;
        o.ki = 0.0;
        o.kd = 0.0;
        o.lo = -d.w_max;
        o.hi = d.w_max;

        PIDParams &n = inner[i].params;
        n.kp = gains.k_w(i);
        n.ki = gains.k_wi(i);
        n.kd = 0.0; // the rate loop has no D, the outer loop's damping is the rate feedback
        n.lo = -tau_max;
        n.hi = tau_max;
        n.anti_windup = AntiWindup::Clamp;
    }
}

void AttitudeController::step(double t, const AttitudeCommand &cmd, const Vector4d &q,
                              const Vector3d &w, SpacecraftControls &u) {
    if (!clock.ready(t)) {
        return;
    }
    const double dt = clock.Ts;

    const Vector4d q_e = cmd.mode == AttitudeCommand::Mode::Point
                             ? pointing_error_quat(q, cmd.b, cmd.d_i)
                             : attitude_error_quat(cmd.q_d, q);

    e = 2.0 * sign * q_e.tail<3>();

    for (int i = 0; i < 3; ++i) {
        // r = 0 and y = e makes the outer pid give -k_theta e
        w_c(i) = outer[i].update(0.0, e(i), dt);
        tau(i) = inner[i].update(w_c(i), w(i), dt);
    }

    // the body feels -A tau_m, so to get tau on the body
    Vector3d tau_m = -A_inv * tau;

    // then it needs to be clamped, so this just scales it down by the max
    const double worst = tau_m.cwiseAbs().maxCoeff();
    if (worst > tau_max && worst > 0.0) {
        tau_m *= tau_max / worst;
    }
    u.wheel_torque = tau_m;
}

void AttitudeController::reset() {
    for (int i = 0; i < 3; ++i) {
        outer[i].reset();
        inner[i].reset();
    }
    sign = 1.0;
    e.setZero();
    w_c.setZero();
    tau.setZero();
}

} // namespace BSIM
