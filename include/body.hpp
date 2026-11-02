#pragma once
#include <Eigen/Dense>
#include <memory>

#include "air_system.hpp"
#include "body_state.hpp"
#include "force_models.hpp"
#include "force_models/fixed_wing_model.hpp"
#include "sensors.hpp"

namespace uav {
using Eigen::Ref;
using Eigen::Vector3d;
using Eigen::Vector4d;

// --- body ---
// a simulated body
// everything is in SI
struct Body {
    double m = 0.0;                // mass
    Vector4d I = Vector4d::Zero(); // moment of inertia (the four components left after the
                                   // symmetry assumptions, Ix, Iy, Iz, Ixz)

    BodyState state; // carries the full 13 dim state

    // the turbulence state for the body (pointer so any Gust model can go there)
    // (null means no turbulence)
    std::unique_ptr<Gust> gust = std::make_unique<DrydenGust>();
    AirState air; // the air around the body this step, it's set by the World before the bodies step

    Sensors sensors; // onboard sensors readings, gets sampled by the World after the bodies step

    // this gets overwritten by subclasses to implement different craft and their specific forces
    virtual Wrench wrench(const BodyState &x, const AirState &air) const = 0;

    // the gust pointer makes bodies move only, and declaring the destructor would otherwise stop
    // the moves being generated
    Body() = default;
    Body(Body &&) = default;
    Body &operator=(Body &&) = default;
    virtual ~Body() = default;
};

// --- body implementations ---
// a basic drone, for now just gravity plus the forces and moments from the keyboard
// uses the zagi's mass properties for now (Table E.1 in appendix E)
struct Drone : Body {
    Wrench input; // the forces and moments the keyboard is asking for (body frame)

    Drone();

    // reads the keyboard into input, call once per frame before stepping
    void read_keyboard();

    Vector3d keyboard_force() const;
    Vector3d keyboard_moment() const;

    // ignores the air, it's only gravity and the keyboard
    Wrench wrench(const BodyState &x, const AirState &) const override;
};

// a fixed wing aircraft with the full chapter 4 force model
// (gravity, longitudinal and lateral aero, propulsion), flying through air defined by AirState
struct FixedWing : Body {
    FixedWingParams params;
    FixedWingControls controls; // control derivatives
    FixedWingControls trim;     // static control values for trim (can equally be 0)

    FixedWing(double m, const Vector4d &I, const FixedWingParams &params);

    // w/s elevator (nose down/up), a/d aileron (roll left/right), q/e rudder (yaw left/right),
    // up/down arrows throttle
    void read_keyboard(double dt);

    Wrench wrench(const BodyState &x, const AirState &air) const override;
};

// --- factory methods ---
// generates a fixed wing with the same characteristics as the zagi flying wing
// everything is from Table E.1 in appendix E
FixedWing generate_zagi_flying_wing();

// generates an imaginary sport plane (not from the book) with the same force model but much more
// control authority, strong yaw stability, very little adverse yaw and about 2x thrust to weight
FixedWing generate_sport_drone();
} // namespace uav
