#include "body_system.hpp"
#include "dynamics.hpp"

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
    return aircraft_derivatives(b.m, b.J, x, F, T);
}

// --- spacecraft ---
SpacecraftBodySystem::Derivative SpacecraftBodySystem::derivatives(const Spacecraft &s,
                                                                   const SpacecraftState &x) const {
    return spacecraft_derivatives(s.params, s.env, s.controls, x);
}

} // namespace BSIM
