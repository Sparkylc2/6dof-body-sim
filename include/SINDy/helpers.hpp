#pragma once

#include "body_system.hpp"
#include "frames.hpp"
#include <Eigen/Dense>

namespace BSIM::SINDy {
using Eigen::MatrixXd;
using Eigen::Vector4d;

// returns [u, w, q, theta]
inline Vector4d to_longitudinal(const BSIM::AircraftState &s) {
    const double theta = quat_to_euler(s.quat()).y();
    return {s.state[BSIM::IDX::VEL],       // u
            s.state[BSIM::IDX::VEL + 2],   // w
            s.state[BSIM::IDX::RATES + 1], // q
            theta};
}

// vel = [u, 0, w], rates = (0, q, 0), pos = (0, 0, -100), quat (0, theta, 0)
inline AircraftState from_longitudinal(const Vector4d &x) {
    BSIM::AircraftState state;
    // vel = (u, 0, w)
    state.vel() = Vector3d(x.x(), 0, x.y());
    state.rates() = Vector3d(0, x.z(), 0);
    state.quat() = euler_to_quat(0, x.w(), 0);
    state.pos() = Vector3d(0, 0, -100);
    return state;
};

// the actual [u_dot, w_dot, q_dot, theta_dot] at the current state
inline Vector4d longitudinal_derivative(const FixedWing &b) {
    const UAVBodySystem sys;
    const auto d = sys.derivatives(b, b.state);

    const double phi = quat_to_euler(b.state.quat()).x();
    const double q = b.state.rates().y();
    const double r = b.state.rates().z();

    // theta_dot from eqn 3.16's middle row, which is just q when wings level
    return {d(IDX::VEL), d(IDX::VEL + 2), d(IDX::RATES + 1), q * std::cos(phi) - r * std::sin(phi)};
}

struct Diffed {
    MatrixXd X, Xdot;
};

// central difference on a trajectory (first and last idx removed since theres no neighbour)
inline Diffed central_diff(const MatrixXd &X, double dt) {
    Diffed out;
    const int m = X.rows();

    out.X = X.middleRows(1, m - 2); // everything but first and last

    // row k of Xdot is x[k+2] - X[k]
    out.Xdot = (X.bottomRows(m - 2) - X.topRows(m - 2)) / (2.0 * dt);
    return out;
}
} // namespace BSIM::SINDy
