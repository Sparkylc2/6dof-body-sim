#pragma once
#include "body_state.hpp"
#include "force_models/spacecraft_model.hpp"

namespace BSIM {

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

} // namespace BSIM
