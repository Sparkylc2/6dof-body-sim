#pragma once
#include "air_system.hpp"
#include "body_state.hpp"
#include "control.hpp"
#include "force_models/fixed_wing_model.hpp"

namespace uav {

// what the autopilot gets told to do
// this is just a sketch from my year 2 lab
// (we did controllers for airspeed, alt, and pitch)
struct AutopilotCommands {
    double Va_c = 0.0;  // airspeed
    double h_c = 0.0;   // altitude
    double phi_c = 0.0; // pitch
};

// what the autopilot can look at
// could be read from both true state and from an estimate built off of the sensors
struct ControlState {
    double pn = 0.0, pe = 0.0, h = 0.0;       // position, h is up
    double Va = 0.0, alpha = 0.0, beta = 0.0; // air data
    double phi = 0.0, theta = 0.0, psi = 0.0; // attitude
    double chi = 0.0, Vg = 0.0;               // course and ground speed
    double p = 0.0, q = 0.0, r = 0.0;         // body rates
    // possible: use accelerometers? later for g-limiting
};

// true state of the body
ControlState true_control_state(const BodyState &x, const AirState &air);

// --- autopilot ---
// a discrete controller that turns the reference autopilot commands and state into control inputs
// specific autopilots will fill in update
struct Autopilot {
    SampleClock clock;          // how often it runs, clock.Ts is the dt for the loops inside
    FixedWingControls controls; // the last control output, stays between samples

    // samples if any are ready at t, otherwise just keeps the last control outputs
    const FixedWingControls &step(double t, const AutopilotCommands &cmd, const ControlState &x);

    // one sample of the controller Ts time after the last sample
    virtual FixedWingControls update(const AutopilotCommands &cmd, const ControlState &x) = 0;

    // just puts every loop back to its starting state
    virtual void reset() = 0;

    virtual ~Autopilot() = default;
};

} // namespace uav
