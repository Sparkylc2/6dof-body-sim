#include "body_system.hpp"
#include "frames.hpp"

namespace uav {

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
    const Vector4d I = b.I;

    // factors in the final expression (B.4) for body angular accel
    const double g = 1.0 / (I.x() * I.z() - I.w() * I.w());
    const double g1 = g * I.w() * (I.x() - I.y() + I.z());
    const double g2 = g * (I.z() * I.z() - I.z() * I.y() + I.w() * I.w());
    const double g3 = g * I.z();
    const double g4 = g * I.w();
    const double g5 = (I.z() - I.x()) / I.y();
    const double g6 = I.w() / I.y();
    const double g7 = g * (I.x() * I.x() - I.y() * I.x() + I.w() * I.w());
    const double g8 = g * I.x();

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
    ang_accel_b.y() += T.y() / I.y();
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

} // namespace uav
