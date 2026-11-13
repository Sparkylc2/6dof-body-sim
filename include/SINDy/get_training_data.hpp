#pragma once

#include "SINDy/helpers.hpp"
#include "body_system.hpp"
#include "trim.hpp"
#include "world.hpp"
#include <Eigen/Dense>

namespace BSIM::SINDy {
using Eigen::MatrixXd;
using Eigen::Vector4d;
using Eigen::VectorXd;

template <int N> struct Trajectory {
    VectorXd t;         // m measurements
    MatrixXd X;         // m x 4 (4 state variables)
    MatrixXd Xdot;      // m x 4 (from central difference)
    MatrixXd Xdot_true; // m x 4 (from exact derivatives)
    void resize(int m) {
        t.resize(m);
        X.resize(m, N);
        Xdot.resize(m, N);
        Xdot_true.resize(m, N);
    }
};

// makes a new plane, sets control to trim, disables turbulence
// gives a little perturbation, steps the World and records
//
// starts at x0 + dx0 (state [u, w, q, theta])
inline Trajectory<4> simulate_longitudinal(const LongitudinalTrim &trim, const Eigen::Vector4d &dx0,
                                           double T, double dt) {

    // initial state (from trim)
    const Vector4d x0 = to_longitudinal(trim.state());

    // world + plane init
    BSIM::World world;
    auto plane = std::make_unique<BSIM::FixedWing>(generate_zagi_flying_wing());
    BSIM::FixedWing &b = *plane;

    // world config
    world.bodies.bodies.push_back(std::move(plane));
    world.atmosphere.turbulence = false;

    // body config
    b.state = from_longitudinal(x0 + dx0);
    b.trim = trim.controls;
    b.controls = trim.controls;

    // trajectory
    const int m = static_cast<int>(std::lround(T / dt));
    Trajectory<4> traj;
    traj.resize(m);

    for (int i = 0; i < m; i++) {

        // collect the measured data
        traj.t(i) = i * dt;
        traj.X.row(i) = (to_longitudinal(b.state) - x0).transpose();
        traj.Xdot_true.row(i) = longitudinal_derivative(b).transpose();

        world.step(dt);
    };

    // doing central diff removes first and last indices
    Diffed diff = central_diff(traj.X, dt);
    traj.t = traj.t.middleRows(1, m - 2).eval();
    traj.Xdot_true = traj.Xdot_true.middleRows(1, m - 2).eval();
    traj.Xdot = diff.Xdot;
    traj.X = diff.X;

    return traj;
};

} // namespace BSIM::SINDy
