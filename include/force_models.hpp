#pragma once
#include "air_system.hpp"
#include <Eigen/Dense>
namespace uav {
using Eigen::Vector3d;
using Eigen::Vector4d;

// gravitational acceleration used for everything (forces and sensors)
constexpr double GRAVITY = 9.81;

// a force and moment acting on a body, both in the body frame
struct Wrench {
    Vector3d F = Vector3d::Zero(); // force [fx, fy, fz]
    Vector3d M = Vector3d::Zero(); // moment [l, m, n]

    Wrench &operator+=(const Wrench &o) {
        F += o.F;
        M += o.M;
        return *this;
    }
};

// airspeed, angle of attack, sideslip and the density of the air the body is in
struct AirData {
    double Va = 0.0;
    double alpha = 0.0;
    double beta = 0.0;
    double rho = RHO_DEFAULT;
};

// air data from the body ground velocity [u, v, w], the attitude and the air around the body
// the wind (steady rotated into body, plus gusts) is taken off to get [u_r, v_r, w_r]
// from eqn 2.8
AirData air_data(const Vector3d &vel_ground_b, const Vector4d &quat, const AirState &air);

} // namespace uav
