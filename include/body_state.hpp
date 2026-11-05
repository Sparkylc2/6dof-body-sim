#pragma once
#include <Eigen/Dense>

namespace BSIM {
using Eigen::Ref;
using Vector1d = Eigen::Vector<double, 1>;
using Eigen::Vector3d;
using Eigen::Vector4d;

template <int N> using StateVector = Eigen::Matrix<double, N, 1>;

// the layout every body shares, anything specific to a type of body goes after the rates
enum IDX {
    POS = 0,   // position
    VEL = 3,   // velocity
    QUAT = 6,  // attitude
    RATES = 10 // rotation rate
};

enum SPACECRAFT_IDX {
    S_REACTION_WHEEL_MOM = 13, // reaction wheel momentum [h_wx, h_wy, h_wz]
    S_MASS = 16                // mass [m]
};

// --- rigid body state ---
// the part of the state body has, [pos | vel | quat | rates | extras...]
// specific bodies inherit this with their own size, and add accessors for whatever extra state
// they put after the rates
template <int N> struct RigidBodyState {
    static_assert(N >= 13, "need at least pos, vel, quat and rates");
    static constexpr int SIZE = N;

    StateVector<N> state = StateVector<N>::Zero();

    RigidBodyState() { quat() << 1.0, 0.0, 0.0, 0.0; }

    // returns a mutable view into the state vector
    // position (inertial FoR) [x, y, z]
    Ref<Vector3d> pos() { return state.template segment<3>(IDX::POS); }
    // velocity, which frame it's in depends on the body
    Ref<Vector3d> vel() { return state.template segment<3>(IDX::VEL); }
    // attitude of body relative to the inertial frame [e0, e1, e2, e3]
    Ref<Vector4d> quat() { return state.template segment<4>(IDX::QUAT); }
    // angular velocity (body FoR) [p, q, r]
    Ref<Vector3d> rates() { return state.template segment<3>(IDX::RATES); }

    // const versions
    Ref<const Vector3d> pos() const { return state.template segment<3>(IDX::POS); }
    Ref<const Vector3d> vel() const { return state.template segment<3>(IDX::VEL); }
    Ref<const Vector4d> quat() const { return state.template segment<4>(IDX::QUAT); }
    Ref<const Vector3d> rates() const { return state.template segment<3>(IDX::RATES); }
};

// --- state implementations ---
// aircraft, layout is the same as B.6
// [pn pe pd | u v w | e0 e1 e2 e3 | p q r]
// (vel is the ground velocity in the body frame)
struct AircraftState : RigidBodyState<13> {};

// spacecraft, the wheels and mass go on the end
// [x y z | vx vy vz | e0 e1 e2 e3 | p q r | h_wx h_wy h_wz | m]
// (pos and vel are both in the inertial frame)
struct SpacecraftState : RigidBodyState<17> {
    // each reaction wheel's total momentum along its own spin axis [h_1, h_2, h_3],
    // i'm probably going to keep the wheels aligned with the body axes so its just the same
    // as the body components
    Ref<Vector3d> wheel_mom() { return state.segment<3>(S_REACTION_WHEEL_MOM); }
    // total mass
    double &mass() { return state(S_MASS); }

    Ref<const Vector3d> wheel_mom() const { return state.segment<3>(S_REACTION_WHEEL_MOM); }
    double mass() const { return state(S_MASS); }
};

} // namespace BSIM
