#pragma once
#include <Eigen/Dense>

#include "air_system.hpp"
#include "body.hpp"
#include "body_state.hpp"
#include "force_models/fixed_wing_model.hpp"

namespace BSIM {
using Eigen::Vector3d;

// --- longitudinal trim (chapter 5) ---
// wings level, straight flight at a constant airspeed and flight path angle in still air
// (gamma = 0 is level flight, positive gamma climbs)
struct LongitudinalTrim {
    double Va = 0.0;            // airspeed
    double gamma = 0.0;         // flight path angle
    double alpha = 0.0;         // theta = alpha + gamma
    FixedWingControls controls; // de and dt get solved

    double residual = 0.0; // ||[u_dot, w_dot, q_dot]|| at the solution
    bool converged = false;

    double theta() const { return alpha + gamma; }

    // trimmed state
    // [u, v, w] = Va [cos(alpha), 0, sin(alpha)] and euler angles (0, theta, 0)
    AircraftState state(const Vector3d &pos = Vector3d::Zero()) const;
};

// finds [alpha, de, dt] so that u_dot = w_dot = q_dot = 0 with theta = alpha + gamma
// uses newton's method
LongitudinalTrim trim_longitudinal(const FixedWing &plane, double Va, double gamma = 0.0,
                                   double rho = RHO_DEFAULT);

} // namespace BSIM
