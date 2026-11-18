#pragma once
#include <Eigen/Dense>
#include <memory>
#include <vector>

#include "actuators.hpp"
#include "air_system.hpp"
#include "body_state.hpp"
#include "force_models.hpp"
#include "force_models/fixed_wing_model.hpp"
#include "force_models/rocket_model.hpp"
#include "force_models/spacecraft_model.hpp"
#include "sensors.hpp"

// --- craft ---
// every vehicle the sim flies: the structs holding a craft's parameters, controls, actuators and
// state together, plus the factories for the built in ones
// (the physics for each is in force_models/ and dynamics.hpp, these just hold the pieces)
namespace BSIM {

using Eigen::Ref;
using Eigen::Vector3d;
using Eigen::Vector4d;

// --- body ---
// a simulated body
// everything is in SI
struct Body {
    double m = 0.0;                // mass
    Vector4d J = Vector4d::Zero(); // moment of inertia (the four components left after the
                                   // symmetry assumptions, Jx, Jy, Jz, Jxz)

    using State = AircraftState; // which state the body carries (what the BodySystem steps)
    State state;                 // carries the full 13 dim state

    // the turbulence state for the body (pointer so any Gust model can go there)
    // (null means no turbulence)
    std::unique_ptr<Gust> gust = std::make_unique<DrydenGust>();
    AirState air; // the air around the body this step, it's set by the World before the bodies step

    Sensors sensors; // onboard sensors readings, gets sampled by the World after the bodies step

    // this gets overwritten by subclasses to implement different craft and their specific forces
    virtual Wrench wrench(const AircraftState &x, const AirState &air) const = 0;

    // moves any actuators towards their commands, the World calls it before the physics step
    // (nothing by default, so bodies without actuator models act on their inputs instantly)
    virtual void step_actuators(double /*dt*/) {}

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
    // the forces and moments being asked for (body frame)
    Wrench input;

    Drone();

    // ignores the air, it's only gravity and the keyboard
    Wrench wrench(const AircraftState &x, const AirState &) const override;
};

// a fixed wing aircraft with the full chapter 4 force model
// (gravity, longitudinal and lateral aero, propulsion), flying through air defined by AirState
struct FixedWing : Body {
    FixedWingParams params;
    FixedWingControls controls; // control derivatives
    FixedWingControls trim;     // static control values for trim (can equally be 0)

    // the aileron servo, controls.da is then the command and aileron.pos is what the aero sees
    // (ideal by default, which is the same as having no actuator)
    Actuator aileron;

    FixedWing(double m, const Vector4d &J, const FixedWingParams &params);

    // uses the current controls (with the aileron where its actuator has got to)
    Wrench wrench(const AircraftState &x, const AirState &air) const override;

    void step_actuators(double dt) override;
    // same force model but with any controls, so things like the trim solver can try out
    // different deflections without touching the body
    Wrench wrench(const AircraftState &x, const AirState &air, const FixedWingControls &u) const;
};

// --- factory methods ---
// generates a fixed wing with the same characteristics as the zagi flying wing
// everything is from Table E.1 in appendix E
FixedWing generate_zagi_flying_wing();

// generates an imaginary sport plane (not from the book) with the same force model but much more
// control authority, strong yaw stability, very little adverse yaw and about 2x thrust to weight
FixedWing generate_sport_drone();

// --- spacecraft ---
// a rigid spacecraft with reaction wheels, thrusters and a fuel tank, in orbit around the earth
// everything is in SI
struct Spacecraft {
    using State = SpacecraftState; // which state it carries (what the BodySystem steps)

    SpacecraftParams params;
    SpacecraftEnvironment env;   // which environment effects act on it
    SpacecraftControls controls; // what the actuators are being asked to do

    State state; // [r | v | q | w | h | m], r and v are inertial (earth centred)
};

// --- factory methods ---
// generates a made up 6U cubesat (about 10 x 23 x 37 cm, 12 kg with 1 kg of fuel) with three
// reaction wheels on the body axes and one main thruster pushing along +z through the cm,
// starting in a 500 km circular orbit at 51.6 deg inclination, not rotating and full of fuel
Spacecraft generate_cubesat();

// --- rocket ---
// a variable mass rocket with a gimballed engine and grid fins
// the physics is all in rocket_derivatives (dynamics.hpp), this just holds the vehicle together
// with its actuators
struct Rocket {
    using State = RocketState;

    RocketParams params;
    RocketEnvironment env;

    // what the controller asks for (the gimbal and fin angles here are commands, the physics sees
    // where the actuators have actually got to, see controls())
    RocketControls commands;

    Actuator gimbal_p, gimbal_y; // pitch and yaw gimbal servos
    std::vector<Actuator> fins;  // one servo per grid fin

    State state; // [p | v | q | w | m], p and v are the cm's, NED

    // moves the gimbal and fins towards their commands over dt (call before the physics step)
    // (the fin servos are held at 0 while the fins are stowed)
    void step_actuators(double dt);

    // the commands with the actuator positions swapped in (what goes into rocket_derivatives)
    // (an ideal actuator is the same as having none, so its command goes straight through)
    RocketControls controls() const;
};

// --- factory methods ---
// a stand in for ICLR's Ember (7.5 m x 0.4 m, see data/openrocket): lox / ethanol with a
// throttleable gimballed engine, four grid fins a bit above the middle (stowed until
// commands.grid_fins_deployed is set), OpenRocket's drag curves, sitting on the pad nose up and
// full of propellant. the masses and engine are placeholders until the real numbers are in
Rocket generate_example_rocket();

} // namespace BSIM
