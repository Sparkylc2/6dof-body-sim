#pragma once
#include "force_models.hpp"
#include <Eigen/Dense>
#include <utility>
#include <vector>

namespace BSIM {
using Eigen::Matrix3d;
using Eigen::Vector3d;

// --- earth ---
constexpr double MU_EARTH = 3.986004418e14; // GM (m^3/s^2)
constexpr double R_EARTH = 6.378137e6;      // equatorial radius (m)
constexpr double J2_EARTH = 1.08262668e-3;  // oblateness
constexpr double G0 = 9.80665;              // standard gravity (just for Isp)

// --- spacecraft parameters ---
// pushes the spacecraft along dir (exhaust goes the other way)
struct Thruster {
    Vector3d pos = Vector3d::Zero();  // where it's mounted (in body FoR)
    Vector3d dir = Vector3d::UnitZ(); // unit direction it pushes the spacecraft (body)
    double T_max = 0.0;               // full throttle thrust (N)
    double Isp = 0.0;                 // specific impulse (s)
};

// 3 reaction wheels, one per col of A (i might make this 4 later since that seems to be standard)
struct ReactionWheels {
    Matrix3d A = Matrix3d::Identity(); // each col is wheel's spin axis in body FoR
    double J_w = 0.0;     // wheel inertia about its spin axis (assuming identical wheels)
    double tau_max = 0.0; // torque limit
    double h_max = 0.0;   // momentum limit (ie top speed)
};

// the empty spacecraft plus a fuel tank
// the body frame's origin is the empty com, when fuel burns the cm will move around in that frame
struct SpacecraftParams {
    double m_empty = 0.0;
    Matrix3d J_empty = Matrix3d::Zero(); // full inertia tensor about the dry cm, wheels locked

    Vector3d r_tank = Vector3d::Zero(); // fuel tank centre (body FoR)
    double R_tank = 0.0;                // tank radius (i'm just assuming it's a solid sphere)

    ReactionWheels wheels;
    std::vector<Thruster> thrusters;
};

// which parts of the environment are switched on (mostly so things can be tested on their own)
struct SpacecraftEnvironment {
    bool gravity = true;
    bool j2 = true;
    bool gravity_gradient = true;
};

// what the actuators are being asked to do, held for the whole step
struct SpacecraftControls {
    std::vector<double> throttle;             // [0, 1] ctrl input per thruster
    Vector3d wheel_torque = Vector3d::Zero(); // motor torque per wheel, body will feel -A * this
};

// --- mass properties ---
// everything about the mass that changes as fuel burns
struct MassProperties {
    double m = 0.0;
    Vector3d r_cm = Vector3d::Zero(); // where the cm currently is (body)
    Matrix3d J = Matrix3d::Zero();    // inertia about the current cm, wheels locked
    Matrix3d J_s = Matrix3d::Zero();  // same with the wheels' spin axis inertia taken out
};

// mass, cm and inertia for the current total mass (fuel = m - m_dry)
MassProperties mass_properties(const SpacecraftParams &p, double m);

// --- gravity ---
// point mass earth, -mu r / |r|^3
Vector3d two_body_gravity(const Vector3d &r);

// the extra acceleration from earth's equatorial bulge
Vector3d j2_gravity(const Vector3d &r);

// gravitational potential per unit mass (gravity = -grad of this) (so i can check energy)
double gravity_potential(const Vector3d &r, bool j2);

// torque from gravity pulls harder on the near side,
// r_b is the position in body axes and J is about the cm
Vector3d gravity_gradient_torque(const Vector3d &r_b, const Matrix3d &J);

// --- actuators ---
// total thrust force and torque about the cm (body FoR) + the fuel burn rate
struct ThrusterOutput {
    Wrench wrench;
    double m_dot = 0.0;
};
ThrusterOutput thrust(const SpacecraftParams &p, const std::vector<double> &throttle,
                      const Vector3d &r_cm, double m);

// the clamped motor torque after applying torque and momentum limits
Vector3d wheel_motor_torque(const ReactionWheels &w, const Vector3d &tau_cmd, const Vector3d &h);

// each wheel's speed relative to the body from its momentum, L / J_w - a dot w
Vector3d wheel_speeds(const ReactionWheels &w, const Vector3d &h, const Vector3d &omega);

// --- orbits ---
// position and velocity (inertial) for a circular orbit at altitude starting on the x axis
std::pair<Vector3d, Vector3d> circular_orbit(double altitude, double inclination);

} // namespace BSIM
