#pragma once
#include "control/blocks.hpp"
#include "force_models/fixed_wing_model.hpp"
#include "sensors.hpp" // DEG
#include "trim.hpp"

namespace BSIM {

// --- longitudinal autopilot ---
// super basic for now, the gains are hand tuned rather than from the design equations in 6.4
//   altitude hold -> pitch command (pi)
//   pitch hold    -> elevator (p on the error, d on the measured pitch rate)
//   airspeed hold -> throttle (p)
// everything is added on top of the trim values, so with no errors it js flies trim
//
// assumes Cm_de < 0 (+de is nose down), which is true for both the zagi and the sport drone
struct LongitudinalGains {
    double k_h = 0.03;   // altitude error -> pitch [rad/m]
    double ki_h = 0.005; // integral of altitude error -> pitch [rad/(m s)]
    double k_theta = 1.5; // pitch error -> elevator
    double k_q = 0.15;    // pitch rate -> elevator (damping)
    double k_V = 0.1;     // airspeed error -> throttle [1/(m/s)]

    double theta_min = -15.0 * DEG, theta_max = 20.0 * DEG; // pitch command limits
    double int_max = 10.0 * DEG; // how much pitch the integrator can add (anti windup)
    double de_max = 30.0 * DEG;
};

struct LongitudinalAutopilot {
    SampleClock clock; // clock.Ts is the dt for the loops
    LongitudinalGains gains;

    double theta_trim = 0.0;
    FixedWingControls trim; // de and dt get used

    double h_int = 0.0;   // integral of the altitude error
    double theta_c = 0.0; // last pitch cmd

    // the trim the loops sit on top of
    void set_trim(const LongitudinalTrim &t);

    // if a sample is ready at t, updates de and dt in u
    // h_c and Va_c are the commands, the rest is what the autopilot receives
    // (either true state or estimate)
    void step(double t, double h_c, double Va_c, double h, double Va, double theta, double q,
              FixedWingControls &u);

    void reset();
};

} // namespace BSIM
