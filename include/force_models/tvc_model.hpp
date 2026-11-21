#pragma once
#include "body_state.hpp"
#include "force_models/rocket_model.hpp"
#include <Eigen/Dense>
#include <array>

// --- hydraulic TVC actuator, and the swinging nozzle it drives ---
// four angled hydraulic pistons push and pull the nozzle about its gimbal
//
// so the controller asks for gimbal angles, inverse kinematics turns them into four piston lengths,
// a proportional loop opens each piston's servo valve, oil flows into cylinders, piston forces
// swing the nozzle, which has its own inertia and the vehicle's loads on it
//
// everything is measured from the gimbal pivot G in body axes
// and the angles mean exactly what they mean in the rocket model
// n_T = C(delta) x with C = R_y(-delta_p) R_z(-delta_y)
//
// the hydraulics are stiff (an oil spring near 100 Hz, millisecond valve lag), so tvc_step sub
// steps

namespace BSIM {
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Matrix42d = Eigen::Matrix<double, 4, 2>;

struct TvcParams {
    // geometry, pistons in a + layout, piston k at roll angle phi0 + 90 k deg
    double x_a = 0.10;    // airframe pins: station from the pivot [m]
    double r_a = 0.17;    //   and radius [m]
    double x_b = -0.15;   // nozzle pins with the nozzle straight: station from the pivot [m]
    double r_b = 0.08;    //   and radius [m]
    double phi0 = 0.0;    // roll angle of piston 0 [rad]
    double stroke = 0.02; // usable travel each side of the straight length [m]

    // hydraulics
    double A_p = 2.0106e-4;   // piston area [m^2] (16 mm bore)
    double P_s = 100e5;       // supply pressure [Pa]
    double Q_max = 4.0212e-5; // valve flow at full opening, no load [m^3/s] (0.2 m/s * A_p)
    double beta = 0.8e9;      // oil bulk modulus [Pa]
    double V_t = 30e-6;       // oil volume per piston, both sides [m^3]
    double C_tp = 2.5e-12;    // leakage past the piston [m^3/s/Pa]
    double tau_v = 3.183e-3;  // valve spool lag [s] (50 Hz)
    double K_p = 300.0;       // position loop gain, valve opening per metre of length error [1/m]

    // nozzle
    double m_n = 6.0;     // [kg]
    double l_n = 0.20;    // its cm behind the pivot [m]
    double I_na = 0.0108; // about its own axis, about its cm [kg m^2]
    double I_nt = 0.125;  // transverse, about its cm [kg m^2]

    // loads the actuator itself adds (section 07)
    double offset_p = 1e-3; // thrust line miss from the pivot, as a pitch moment arm [m]
    double offset_y = 0.0;  //   and yaw [m]
    double friction = 4.0;  // coulomb gimbal friction [N m]
    double damping = 0.0;   // viscous gimbal friction [N m s/rad]

    // integration
    double h = 20e-6; // internal sub step [s]
};

// the actuator's own params with the nozzle taken from the rocket's engine
TvcParams tvc_params_for(const RocketEngine &e);

// the 12 states [delta (2) | delta_dot (2) | P_L (4) | u_v (4)]
//      delta  gimbal angles [pitch, yaw] [rad]
//      P_L    each piston's load pressure, P_1 - P_2 [Pa] (force = A_p P_L)
//      u_v    each valve's spool opening, -1 to 1
constexpr int TVC_NX = 12;
using TvcState = StateVector<TVC_NX>;

// what it gets from outside every step
struct TvcInputs {
    Vector2d delta_cmd = Vector2d::Zero(); // gimbal angles the controller wants [rad]
    Vector2d load = Vector2d::Zero();      // vehicle load on the nozzle [N m] (rocket_nozzle_load)
    double thrust = 0.0;                   // [N], for the thrust offset moment
};

// --- kinematics ---
struct TvcKinematics {
    Vector4d l = Vector4d::Zero();   // piston lengths [m]
    std::array<Vector3d, 4> u;       // unit vectors along each piston, airframe end -> nozzle end
    std::array<Vector3d, 4> r;       // nozzle pins from the pivot [m]
    Matrix42d J = Matrix42d::Zero(); // l_dot = J delta_dot, and gimbal torque = J^T F
};

// airframe pins and the nozzle pins with the nozzle straight
std::array<Vector3d, 4> tvc_anchors(const TvcParams &p);
std::array<Vector3d, 4> tvc_nozzle_pins(const TvcParams &p);

// piston lengths, directions and the jacobian at gimbal angles delta
TvcKinematics tvc_kinematics(const TvcParams &p, const Vector2d &delta);

// gimbal angles from measured piston lengths, gauss newton from delta0
// (residual is the leftover length misfit [m])
Vector2d tvc_forward_kinematics(const TvcParams &p, const Vector4d &l_meas, Vector2d delta0,
                                int iterations = 5, double *residual = nullptr);

// --- dynamics ---
struct TvcOutputs {
    Vector2d delta_ddot = Vector2d::Zero(); // nozzle angular acceleration [rad/s^2]
    Vector4d l = Vector4d::Zero();          // piston lengths [m]
    Vector4d l_cmd = Vector4d::Zero();      // what each piston is being asked for [m]
    Vector4d F = Vector4d::Zero();          // piston forces, + pushes the nozzle pin outwards [N]
    Vector4d Q = Vector4d::Zero();          // valve flows [m^3/s]
};

// the nozzle's 2x2 mass matrix E^T I_G E about its gimbal axes
Eigen::Matrix2d tvc_mass_matrix(const TvcParams &p, const Vector2d &delta);

// the 12 state derivatives
TvcState tvc_derivatives(const TvcParams &p, const TvcState &s, const TvcInputs &in,
                         TvcOutputs *out = nullptr);

// moves the actuator on by dt (rk4 at p.h inside)
// out (optional) is at the end of the step, with delta_ddot the average over it,
// (delta_dot_end - delta_dot_start) / dt
TvcState tvc_step(const TvcParams &p, const TvcState &s, const TvcInputs &in, double dt,
                  TvcOutputs *out = nullptr);

// the state holding still at delta with no load
TvcState tvc_rest(const Vector2d &delta = Vector2d::Zero());

// --- the swinging nozzle's effect on the vehicle ---
//      tail wags dog       -B delta_ddot, the reaction to accelerating N
//      coriolis pair       -w x B delta_dot - (dJ/ddelta delta_dot) w (the two halves go together)
//      changing geometry   -B_dot delta_dot, and the part of J_dot from N swinging
//      swinging exit       -m_dot_e rho_e x rho_e_dot, the nozzle throwing the exhaust sideways
// the cm is kept on the body axis

// N's orientation relative to the core, C_bn = R_y(-delta_p) R_z(-delta_y), so n_T = C_bn x
Eigen::Matrix3d gimbal_C_bn(double delta_p, double delta_y);

// N's angular velocity relative to the core is E(delta) delta_dot (outer ring then inner ring)
Matrix32d gimbal_E(double delta_p);

// N's inertia about the pivot G, in body axes, with N at delta
Eigen::Matrix3d gimbal_I_G(const RocketEngine &e, double delta_p, double delta_y);

// N's inertia about the system cm (at station x_cm, on the axis), with N at delta
Eigen::Matrix3d gimbal_J_N(const RocketEngine &e, double x_cm, double delta_p, double delta_y);

// H_rel = B delta_dot, the angular momentum from N swinging, B = (I_G - m_n [a]x [b]x) E
Matrix32d gimbal_B(const RocketEngine &e, double x_cm, double delta_p, double delta_y);

// every swinging nozzle moment about the cm (body frame), GimbalDynamics is in rocket_model.hpp
//   x_cm_dot is how fast the cm is sliding along the axis as propellant burns
//   rho_e is the nozzle exit from the cm (EngineOutput::rho_e)
GimbalDynamics rocket_gimbal_dynamics(const RocketParams &p, const RocketControls &u,
                                      const Vector3d &w, double x_cm, double x_cm_dot,
                                      double m_dot_e, const Vector3d &rho_e);

// torque the gimbal actuators have to supply [pitch, yaw], positive the same
// way as the gimbal angles.
// needs the body's angular acceleration and specific force (non gravitational force / m, body), so
// it comes after the rates are solved for
Vector2d rocket_hinge_torque(const RocketParams &p, const RocketControls &u, const Vector3d &w,
                             const Vector3d &w_dot, double x_cm, const Vector3d &f_body);

// the load the vehicle puts on N about its gimbal axes [pitch, yaw]
// (the -B^T w_dot, pendulum and velocity terms of the nozzle row), positive pushes the gimbal
// angles up.
Vector2d rocket_nozzle_load(const RocketParams &p, const RocketControls &u, const Vector3d &w,
                            const Vector3d &w_dot, double x_cm, const Vector3d &f_body);

} // namespace BSIM
