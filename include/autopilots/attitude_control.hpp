#pragma once
#include "control.hpp"
#include "force_models/spacecraft_model.hpp"
#include <array>

namespace BSIM {
using Eigen::Matrix3d;
using Eigen::Vector3d;
using Eigen::Vector4d;

// basically anything on here has derivations on page 10 of my attached working out

// --- attitude error ---
// the error quaternion is the body attitude measured from the desired one,
// q_e = q_d (x) q = [cos(phi/2), sin(phi/2) n] with n in body axes
Vector4d attitude_error_quat(const Vector4d &q_d, const Vector4d &q);

// same but for pointing a body vector b along the inertial direction d_i (so roll about b is free)
// q_d is built from q by the shortest rotation taking b onto d, making q_e = q_s*q
//   q_s = [1 + b dot d, b cross d] / sqrt(2 (1 + b dot d)) with d = R(q)^T d_i
// q_e(0) >= 0 always
Vector4d pointing_error_quat(const Vector4d &q, const Vector3d &b, const Vector3d &d_i);

// --- commands ---
// either hold a full attitude q_d or just point body axis b at the inertial direction d_i
struct AttitudeCommand {
    enum class Mode { Attitude, Point };
    Mode mode = Mode::Attitude;

    Vector4d q_d = Vector4d(1.0, 0.0, 0.0, 0.0); // desired attitude (Attitude)

    Vector3d b = Vector3d::UnitZ();   // body axis to point (Point)
    Vector3d d_i = Vector3d::UnitX(); // inertial direction to point it at (Point)
};

// --- tuning parameters ---
struct AttitudeDesignParams {
    double w_rate = 1.0;   // inner rate loop bandwidth (rad/s)
    double zeta = 1.0;     // damping of the combined loop, which sets k_theta = w_rate / (4 zeta^2)
    double ki_ratio = 0.1; // inner PI zero at ki_ratio * w_rate (ie a decade below by default)
    double w_max = 0.05;   // slew rate cap (rad/s), keeps it under l_max / J_ii with some margin
};

// --- gains ---
struct AttitudeGains {
    Vector3d k_w = Vector3d::Zero();  // rate P per axis, J_ii w_rate
    Vector3d k_wi = Vector3d::Zero(); // rate I per axis, ki_ratio w_rate k_w
    double k_theta = 0.0;             // attitude P, w_rate / (4 zeta^2)
    double wn = 0.0;                  // nat freq of the combined loop, sqrt(w_rate k_theta)
};

// only the diagonal of J_s is used (i'm assuming this is a cubesat for now)
AttitudeGains attitude_gains(const Matrix3d &J_s, const AttitudeDesignParams &d);

// --- attitude controller ---
// cascade, per body axis:
//   outer -> w_c = sat(-k_theta e, +-w_max)
//   inner -> tau = sat(k_w (w_c - w) + k_wi/s (w_c - w), +-tau_max)
//
// with e = 2 sgn(e_0e) eps_e from the error quaternion, then the wheels get tau_m = -A^-1 tau
// since the body feels -A tau_m
struct AttitudeController {
    SampleClock clock; // clock.Ts is the dt for the loops

    std::array<PID, 3> outer; // attitude P, one per body axis
    std::array<PID, 3> inner; // rate PI, one per body axis

    AttitudeGains gains;
    Matrix3d A_inv = Matrix3d::Identity(); // wheel axes to body, inverted
    double tau_max = 0.0;

    double sign = 1.0;

    // last values (for debugging/whatever, usually good practice xD)
    Vector3d e = Vector3d::Zero();   // attitude error (rad, roughly, about each body axis)
    Vector3d w_c = Vector3d::Zero(); // rate command
    Vector3d tau = Vector3d::Zero(); // body torque asked for

    // sets the gains from the spacecraft's inertia at mass m and the wheel limits
    void design(const SpacecraftParams &p, double m, const AttitudeDesignParams &d);

    // if a sample is ready at t, updates the wheel torque in u
    // q and w are what the controller receives (true state or estimate)
    void step(double t, const AttitudeCommand &cmd, const Vector4d &q, const Vector3d &w,
              SpacecraftControls &u);

    void reset();
};

} // namespace BSIM
