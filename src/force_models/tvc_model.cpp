#include "force_models/tvc_model.hpp"
#include <algorithm>
#include <cmath>

namespace BSIM {

namespace {

constexpr double PI = 3.14159265358979323846;

// the rocket engine with this nozzle, for the shared gimbal helpers
RocketEngine nozzle_of(const TvcParams &p) {
    RocketEngine e;
    e.m_n = p.m_n;
    e.l_n = p.l_n;
    e.I_na = p.I_na;
    e.I_nt = p.I_nt;
    return e;
}

double sgn(double v) { return (v > 0.0) - (v < 0.0); }

// smoothed sign for coulomb friction
constexpr double FRICTION_RATE = 0.01; // [rad/s]

} // namespace

TvcParams tvc_params_for(const RocketEngine &e) {
    TvcParams p;
    p.m_n = e.m_n;
    p.l_n = e.l_n;
    p.I_na = e.I_na;
    p.I_nt = e.I_nt;
    return p;
}

std::array<Vector3d, 4> tvc_anchors(const TvcParams &p) {
    std::array<Vector3d, 4> a;
    for (int k = 0; k < 4; k++) {
        const double phi = p.phi0 + 0.5 * PI * k;
        a[k] = Vector3d(p.x_a, p.r_a * std::cos(phi), p.r_a * std::sin(phi));
    }
    return a;
}

std::array<Vector3d, 4> tvc_nozzle_pins(const TvcParams &p) {
    std::array<Vector3d, 4> b;
    for (int k = 0; k < 4; k++) {
        const double phi = p.phi0 + 0.5 * PI * k;
        b[k] = Vector3d(p.x_b, p.r_b * std::cos(phi), p.r_b * std::sin(phi));
    }
    return b;
}

TvcKinematics tvc_kinematics(const TvcParams &p, const Vector2d &delta) {
    const Eigen::Matrix3d C = gimbal_C_bn(delta.x(), delta.y());
    const Matrix32d E = gimbal_E(delta.x());

    const auto a = tvc_anchors(p);
    const auto b0 = tvc_nozzle_pins(p);

    TvcKinematics k;
    for (int i = 0; i < 4; i++) {
        k.r[i] = C * b0[i];               // pin rotates with the nozzle
        const Vector3d d = k.r[i] - a[i]; // anchor -> pin

        k.l(i) = d.norm();
        k.u[i] = d / k.l(i);
        // l_dot = (r x u) . Omega, Omega = E delta_dot
        k.J.row(i) = k.r[i].cross(k.u[i]).transpose() * E;
    }
    return k;
}

Vector2d tvc_forward_kinematics(const TvcParams &p, const Vector4d &l_meas, Vector2d delta,
                                int iterations, double *residual) {
    for (int it = 0; it < iterations; it++) {
        const TvcKinematics k = tvc_kinematics(p, delta);
        const Vector4d err = l_meas - k.l;

        delta += (k.J.transpose() * k.J).ldlt().solve(k.J.transpose() * err);
    }

    if (residual) {
        *residual = (l_meas - tvc_kinematics(p, delta).l).norm();
    }

    return delta;
}

Eigen::Matrix2d tvc_mass_matrix(const TvcParams &p, const Vector2d &delta) {
    const Matrix32d E = gimbal_E(delta.x());

    return E.transpose() * gimbal_I_G(nozzle_of(p), delta.x(), delta.y()) * E;
}

TvcState tvc_rest(const Vector2d &delta) {
    TvcState s = TvcState::Zero();
    s.segment<2>(0) = delta;
    return s;
}

TvcState tvc_derivatives(const TvcParams &p, const TvcState &s, const TvcInputs &in,
                         TvcOutputs *out) {
    const Vector2d delta = s.segment<2>(0);
    const Vector2d delta_dot = s.segment<2>(2);

    const Vector4d P_L = s.segment<4>(4);
    const Vector4d u_v = s.segment<4>(8);

    const TvcKinematics k = tvc_kinematics(p, delta);
    const Vector4d l_dot = k.J * delta_dot;

    // IK on the cmd, each piston asked for it's length
    const Vector4d l0 = tvc_kinematics(p, Vector2d::Zero()).l;
    Vector4d l_cmd = tvc_kinematics(p, in.delta_cmd).l;

    for (int i = 0; i < 4; i++) {
        l_cmd(i) = std::clamp(l_cmd(i), l0(i) - p.stroke, l0(i) + p.stroke);
    }

    TvcState sd;
    Vector4d Q, F;
    for (int i = 0; i < 4; i++) {
        // valve - merritt's square root law
        const double drop = std::max(1.0 - sgn(u_v(i)) * P_L(i) / p.P_s, 0.0);
        Q(i) = p.Q_max * u_v(i) * std::sqrt(drop);

        // oil
        sd(4 + i) = 4.0 * p.beta / p.V_t * (Q(i) - p.A_p * l_dot(i) - p.C_tp * P_L(i));

        // spool - lags a proportional loop on the length error
        const double u_cmd = std::clamp(p.K_p * (l_cmd(i) - k.l(i)), -1.0, 1.0);
        sd(8 + i) = (u_cmd - u_v(i)) / p.tau_v;

        F(i) = p.A_p * P_L(i);
    }

    // nozzle - M delta_ddot = J^T F + tau_ext
    const Vector2d offset(p.offset_p, p.offset_y);
    Vector2d tau_ext = in.load + in.thrust * offset - p.damping * delta_dot;

    for (int j = 0; j < 2; j++) {
        tau_ext(j) -= p.friction * std::tanh(delta_dot(j) / FRICTION_RATE);
    }

    const Vector2d delta_ddot =
        tvc_mass_matrix(p, delta).ldlt().solve(k.J.transpose() * F + tau_ext);

    sd.segment<2>(0) = delta_dot;
    sd.segment<2>(2) = delta_ddot;

    if (out) {
        out->delta_ddot = delta_ddot;
        out->l = k.l;
        out->l_cmd = l_cmd;
        out->F = F;
        out->Q = Q;
    }

    return sd;
}

TvcState tvc_step(const TvcParams &p, const TvcState &s, const TvcInputs &in, double dt,
                  TvcOutputs *out) {
    const int n = std::max(1, static_cast<int>(std::ceil(dt / p.h)));
    const double h = dt / n;

    TvcState x = s;
    for (int i = 0; i < n; i++) {
        const TvcState k1 = tvc_derivatives(p, x, in);
        const TvcState k2 = tvc_derivatives(p, x + 0.5 * h * k1, in);
        const TvcState k3 = tvc_derivatives(p, x + 0.5 * h * k2, in);
        const TvcState k4 = tvc_derivatives(p, x + h * k3, in);
        x += h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }

    if (out) {
        tvc_derivatives(p, x, in, out);

        if (dt > 0.0) {
            out->delta_ddot = (x.segment<2>(2) - s.segment<2>(2)) / dt;
        }
    }
    return x;
}

// --- the swinging nozzle's effect on the vehicle ---
namespace {

// [v]x, so [v]x u = v x u
Eigen::Matrix3d skew(const Vector3d &v) {
    Eigen::Matrix3d S;
    S << 0.0, -v.z(), v.y(), v.z(), 0.0, -v.x(), -v.y(), v.x(), 0.0;
    return S;
}

// a point mass's inertia about the origin it's measured from, (r.r) I - r r^T
Eigen::Matrix3d point_inertia(double m, const Vector3d &r) {
    return m * (r.squaredNorm() * Eigen::Matrix3d::Identity() - r * r.transpose());
}

// N's inertia about its own cm in its own axes
Eigen::Matrix3d own_inertia(const RocketEngine &e) {
    return Vector3d(e.I_na, e.I_nt, e.I_nt).asDiagonal();
}

// cm -> pivot (a) and pivot -> N's cm (b), body frame
Vector3d pivot_from_cm(const RocketEngine &e, double x_cm) { return {e.x_G - x_cm, 0.0, 0.0}; }
Vector3d n_cm_from_pivot(const RocketEngine &e, double delta_p, double delta_y) {
    return -e.l_n * thrust_direction(delta_p, delta_y);
}

// a change of J or B over a small time step along how things are moving (delta at delta_dot, the cm
// at x_cm_dot), central difference
constexpr double H_STEP = 1e-5; // [s]

} // namespace

Eigen::Matrix3d gimbal_C_bn(double delta_p, double delta_y) {
    const double cp = std::cos(delta_p), sp = std::sin(delta_p);
    const double cy = std::cos(delta_y), sy = std::sin(delta_y);
    Eigen::Matrix3d Ry, Rz; // R_y(-delta_p), R_z(-delta_y)
    Ry << cp, 0.0, -sp, 0.0, 1.0, 0.0, sp, 0.0, cp;
    Rz << cy, sy, 0.0, -sy, cy, 0.0, 0.0, 0.0, 1.0;
    return Ry * Rz;
}

Matrix32d gimbal_E(double delta_p) {
    Matrix32d E;
    E << 0.0, std::sin(delta_p), -1.0, 0.0, 0.0, -std::cos(delta_p);
    return E;
}

Eigen::Matrix3d gimbal_I_G(const RocketEngine &e, double delta_p, double delta_y) {
    const Eigen::Matrix3d C = gimbal_C_bn(delta_p, delta_y);
    return C * own_inertia(e) * C.transpose() +
           point_inertia(e.m_n, n_cm_from_pivot(e, delta_p, delta_y));
}

Eigen::Matrix3d gimbal_J_N(const RocketEngine &e, double x_cm, double delta_p, double delta_y) {
    const Eigen::Matrix3d C = gimbal_C_bn(delta_p, delta_y);
    const Vector3d rho_n = pivot_from_cm(e, x_cm) + n_cm_from_pivot(e, delta_p, delta_y);
    return C * own_inertia(e) * C.transpose() + point_inertia(e.m_n, rho_n);
}

Matrix32d gimbal_B(const RocketEngine &e, double x_cm, double delta_p, double delta_y) {
    const Vector3d a = pivot_from_cm(e, x_cm);
    const Vector3d b = n_cm_from_pivot(e, delta_p, delta_y);
    return (gimbal_I_G(e, delta_p, delta_y) - e.m_n * skew(a) * skew(b)) * gimbal_E(delta_p);
}

GimbalDynamics rocket_gimbal_dynamics(const RocketParams &p, const RocketControls &u,
                                      const Vector3d &w, double x_cm, double x_cm_dot,
                                      double m_dot_e, const Vector3d &rho_e) {
    const RocketEngine &e = p.engine;
    const Vector2d d(u.gimbal_p, u.gimbal_y);
    const Vector2d d_dot(u.gimbal_p_dot, u.gimbal_y_dot);
    const Vector2d d_ddot(u.gimbal_p_ddot, u.gimbal_y_ddot);

    GimbalDynamics g;

    // exit moves at -L_n (Omega_n x n_T) relative to the core (the cm sliding
    // with the propellant is a flow transient)
    const Vector3d Omega_n = gimbal_E(d.x()) * d_dot;
    const Vector3d rho_e_dot = -e.L_n * Omega_n.cross(thrust_direction(d.x(), d.y()));
    g.M_exit = -m_dot_e * rho_e.cross(rho_e_dot);

    if (e.m_n <= 0.0) {
        return g; // massless nozzle so only the exhaust term
    }

    auto dJ_at = [&](double xc, const Vector2d &dd) -> Eigen::Matrix3d {
        return gimbal_J_N(e, xc, dd.x(), dd.y()) - gimbal_J_N(e, xc, 0.0, 0.0);
    };

    g.dJ = dJ_at(x_cm, d);
    g.B = gimbal_B(e, x_cm, d.x(), d.y());
    g.H_rel = g.B * d_dot;
    g.M_twd = -g.B * d_ddot;

    // how dJ and B change over time, split into the part from N swinging (delta moving) and the
    // part from the cm sliding (x_cm moving), central differences along the motion
    const double h = H_STEP;
    const Eigen::Matrix3d dJ_delta =
        (dJ_at(x_cm, d + h * d_dot) - dJ_at(x_cm, d - h * d_dot)) / (2.0 * h);
    const Eigen::Matrix3d dJ_cm =
        (dJ_at(x_cm + h * x_cm_dot, d) - dJ_at(x_cm - h * x_cm_dot, d)) / (2.0 * h);
    const Vector2d dp = d + h * d_dot, dm = d - h * d_dot;
    const Matrix32d B_dot = (gimbal_B(e, x_cm + h * x_cm_dot, dp.x(), dp.y()) -
                             gimbal_B(e, x_cm - h * x_cm_dot, dm.x(), dm.y())) /
                            (2.0 * h);

    g.M_cor = -w.cross(g.H_rel) - dJ_delta * w;
    g.M_geom = -B_dot * d_dot - dJ_cm * w;
    return g;
}

Vector2d rocket_hinge_torque(const RocketParams &p, const RocketControls &u, const Vector3d &w,
                             const Vector3d &w_dot, double x_cm, const Vector3d &f_body) {
    const RocketEngine &e = p.engine;
    if (e.m_n <= 0.0) {
        return Vector2d::Zero();
    }
    const double dp = u.gimbal_p, dy = u.gimbal_y;
    const Vector2d d_dot(u.gimbal_p_dot, u.gimbal_y_dot);
    const Vector2d d_ddot(u.gimbal_p_ddot, u.gimbal_y_ddot);

    const Matrix32d E = gimbal_E(dp);
    const Vector3d Omega_n = E * d_dot;
    const Vector3d Omega_n_dot =
        E * d_ddot + d_dot.x() * d_dot.y() * Vector3d(std::cos(dp), 0.0, std::sin(dp));

    // N's absolute rate and angular acceleration, in core axes
    const Vector3d w_n = w + Omega_n;
    const Vector3d alpha_n = w_dot + Omega_n_dot + w.cross(Omega_n);

    // the pivot's acceleration minus gravity (gravity acts on N's cm, so it cancels like in f)
    const Vector3d a = pivot_from_cm(e, x_cm);
    const Vector3d b = n_cm_from_pivot(e, dp, dy);
    const Vector3d a_G = f_body + w_dot.cross(a) + w.cross(w.cross(a));

    // euler's equation for N about its moving pivot
    const Eigen::Matrix3d I_G = gimbal_I_G(e, dp, dy);
    const Vector3d tau_h = I_G * alpha_n + w_n.cross(I_G * w_n) + e.m_n * b.cross(a_G);

    // each actuator axis takes the part along its own axis: pitch about -y, yaw about -e_y
    const Vector3d e_y(-std::sin(dp), 0.0, std::cos(dp));
    return {-tau_h.y(), -e_y.dot(tau_h)};
}

Vector2d rocket_nozzle_load(const RocketParams &p, const RocketControls &u, const Vector3d &w,
                            const Vector3d &w_dot, double x_cm, const Vector3d &f_body) {
    RocketControls still = u;
    still.gimbal_p_ddot = still.gimbal_y_ddot = 0.0;
    return -rocket_hinge_torque(p, still, w, w_dot, x_cm, f_body);
}

} // namespace BSIM
