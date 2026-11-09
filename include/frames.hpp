// include/uav/frames.hpp
#pragma once
#include <Eigen/Dense>

namespace BSIM {
using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::Vector3d;
using Eigen::Vector4d;

// --- euler ---
// returns the rotation matrix rotating the reference system by euler angles
// based on the content of chapter 2 (2.4 and 2.5 in particular)
Matrix3d rotation_euler(double phi, double theta, double psi);

// rotation matrix bringing body angular velocity to the inertial frame of reference
// from eqn 3.16
Matrix3d body_rates_to_euler_rates(double phi, double theta);

// --- quaternion ---
// returns a rotation matrix that rotates the reference system by the quaternion e
// from eqn B.1
Matrix3d rotation(const Vector4d &e);

// rate of change of the attitude quaternion based on the current attitude and body rates [p, q, r]
// from eqn B.5
Vector4d quat_rate(const Vector4d &e, const Vector3d &w);

// converts euler to quaternion and vice versa
// from eqn B.3
Vector4d euler_to_quat(double phi, double theta, double psi);
Vector3d quat_to_euler(const Vector4d &e);

// hamilton product p (x) q (scalar first)
Vector4d quat_mul(const Vector4d &p, const Vector4d &q);

// conjugate for unit quaternion
Vector4d quat_conj(const Vector4d &q);

} // namespace BSIM
