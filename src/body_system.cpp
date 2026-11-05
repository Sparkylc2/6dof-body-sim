#include "body_system.hpp"
#include "frames.hpp"

namespace BSIM {

// --- generic ---
template <class B> void BodySystem<B>::step(double dt) {
    for (const std::unique_ptr<B> &body : bodies) {
        B &b = *body;

        // slope of the state at x, each system works out its own forces inside derivatives
        auto f = [this, &b](const State &x) { return derivatives(b, x); };

        // --- rk4 ---
        // the state at the start of the step
        const State &x0 = b.state;
        auto offset = [&x0](const Derivative &k, double h) {
            State x = x0;
            x.state += h * k;
            return x;
        };

        // standard rk4
        const Derivative k1 = f(x0);
        const Derivative k2 = f(offset(k1, dt / 2.0));
        const Derivative k3 = f(offset(k2, dt / 2.0));
        const Derivative k4 = f(offset(k3, dt));
        b.state.state += dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

        // the update makes ||e|| larger than 1, so it needs to be normalised
        b.state.quat().normalize();
    }
}

// --- aircraft ---
UAVBodySystem::Derivative UAVBodySystem::derivatives(const Body &b, const AircraftState &x) const {
    // the air is held for the whole step, but the steady wind still gets rotated by each stage's
    // attitude inside the wrench
    const Wrench w = b.wrench(x, b.air);
    return derivatives(b, x, w.F, w.M);
}

UAVBodySystem::Derivative UAVBodySystem::derivatives(const Body &b, const AircraftState &x,
                                                     const Vector3d &F, const Vector3d &T) const {
    // the resulting system derivative
    Derivative x_dot;

    // state (b - body, v - vehicle) frame
    const Vector3d vel_b = x.vel();     // u, v, w
    const Vector3d rates_b = x.rates(); // p, q, r

    // moment of inertia
    const Vector4d J = b.J;

    // factors in the final expression (B.4) for body angular accel
    const double g = 1.0 / (J.x() * J.z() - J.w() * J.w());
    const double g1 = g * J.w() * (J.x() - J.y() + J.z());
    const double g2 = g * (J.z() * J.z() - J.z() * J.y() + J.w() * J.w());
    const double g3 = g * J.z();
    const double g4 = g * J.w();
    const double g5 = (J.z() - J.x()) / J.y();
    const double g6 = J.w() / J.y();
    const double g7 = g * (J.x() * J.x() - J.y() * J.x() + J.w() * J.w());
    const double g8 = g * J.x();

    // --- body-frame translational dynamics ---
    // the -w x v_b term is from differentiating in the rotating frame
    // from eqn B.2
    Vector3d accel_b;
    // acceleration due to rotation of the body
    accel_b.x() = rates_b.z() * vel_b.y() - rates_b.y() * vel_b.z();
    accel_b.y() = rates_b.x() * vel_b.z() - rates_b.z() * vel_b.x();
    accel_b.z() = rates_b.y() * vel_b.x() - rates_b.x() * vel_b.y();

    // body-frame forces
    accel_b += 1 / b.m * F;

    // --- body-frame rotational dynamics ---
    // from eqn B.4
    Vector3d ang_accel_b;

    // gyroscopic terms (i think if thats what they are called)
    ang_accel_b.x() = g1 * rates_b.x() * rates_b.y() - g2 * rates_b.y() * rates_b.z();
    ang_accel_b.y() = g5 * rates_b.x() * rates_b.z() -
                      g6 * (rates_b.x() * rates_b.x() - rates_b.z() * rates_b.z());
    ang_accel_b.z() = g7 * rates_b.x() * rates_b.y() - g1 * rates_b.y() * rates_b.z();

    // body-frame torques
    ang_accel_b.x() += g3 * T.x() + g4 * T.z();
    ang_accel_b.y() += T.y() / J.y();
    ang_accel_b.z() += g4 * T.x() + g8 * T.z();

    // velocity of the body in inertial frame
    x_dot.segment<3>(IDX::POS) = rotation(x.quat()) * x.vel();
    // acceleration of the body in the body frame
    x_dot.segment<3>(IDX::VEL) = accel_b;
    // rate of change of the quaternion going from velocity to inertial frame
    x_dot.segment<4>(IDX::QUAT) = quat_rate(x.quat(), x.rates());
    // angular acceleration of the body in the body frame
    x_dot.segment<3>(IDX::RATES) = ang_accel_b;

    return x_dot;
}

// --- spacecraft ---
SpacecraftBodySystem::Derivative SpacecraftBodySystem::derivatives(const Spacecraft &s,
                                                                   const SpacecraftState &x) const {
    const SpacecraftParams &p = s.params;
    const SpacecraftEnvironment &env = s.env;

    // the resulting system derivative
    Derivative x_dot;

    const Vector3d r = x.pos(); // inertial
    const Vector3d w = x.rates();
    const Vector3d l = x.wheel_mom();
    const Matrix3d R = rotation(x.quat()); // body to inertial
    const MassProperties mp =
        mass_properties(p, x.mass()); // the cm and inertia will move/change with fuel burn

    // thrust force and torque about the current cm, and the fuel burn info
    const ThrusterOutput th = thrust(p, s.controls.throttle, mp.r_cm, x.mass());

    // -- translational accel --- (it's all in the inertial frame so no cross term this time)
    // thrust + gravity + equatorial bulge
    Vector3d accel = R * th.wrench.F / mp.m;
    if (env.gravity)
        accel += two_body_gravity(r);
    if (env.j2)
        accel += j2_gravity(r);

    // --- rotational accel ---
    // thrust + gravity_gradient
    Vector3d tau = th.wrench.M;
    if (env.gravity && env.gravity_gradient)
        tau += gravity_gradient_torque(R.transpose() * r, mp.J);

    // the wheels have a torque about their axis so the satellite will feel
    // -A * tau_m (ie A gets torque from the wheel frame to the body frame)
    const Vector3d tau_m = wheel_motor_torque(p.wheels, s.controls.wheel_torque, l);
    const Matrix3d &A = p.wheels.A;
    const Vector3d l_w = A * l; // wheel momentum in body axes

    // euler with the wheels is:
    // J_s * w_dot + h_w_dot + w cross (J_s*w + h_w) = tau
    // (J would change as fuel is used but i'm assuming it's really small so it's left out)
    const Vector3d w_dot = mp.J_s.inverse() * (tau - A * tau_m - w.cross(mp.J_s * w + l_w));

    x_dot.segment<3>(IDX::POS) = x.vel();
    x_dot.segment<3>(IDX::VEL) = accel;
    x_dot.segment<4>(IDX::QUAT) = quat_rate(x.quat(), w);
    x_dot.segment<3>(IDX::RATES) = w_dot;
    // the motor torque is the change in each wheel's axial momentum by definition so
    x_dot.segment<3>(SPACECRAFT_IDX::S_REACTION_WHEEL_MOM) = tau_m;
    x_dot(SPACECRAFT_IDX::S_MASS) = th.m_dot;

    return x_dot;
}

} // namespace BSIM
