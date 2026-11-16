#include "control/autopilots/longitudinal_autopilot.hpp"
#include <algorithm>

namespace BSIM {

void LongitudinalAutopilot::set_trim(const LongitudinalTrim &t) {
    theta_trim = t.theta();
    trim = t.controls;
}

void LongitudinalAutopilot::step(double t, double h_c, double Va_c, double h, double Va,
                                 double theta, double q, FixedWingControls &u) {
    if (!clock.ready(t)) {
        return;
    }
    const double dt = clock.Ts;
    const LongitudinalGains &g = gains;

    // altitude -> pitch, the integrator is clamped so it can only add +-int_max of pitch
    const double e_h = h_c - h;
    const double int_lim = g.int_max / g.ki_h;
    h_int = std::clamp(h_int + e_h * dt, -int_lim, int_lim);
    theta_c = std::clamp(theta_trim + g.k_h * e_h + g.ki_h * h_int, g.theta_min, g.theta_max);

    // pitch -> elevator, minus on the error since +de pitches the nose down
    // (and + on q so a nose up rate gets pushed back down)
    u.de = std::clamp(trim.de - g.k_theta * (theta_c - theta) + g.k_q * q, -g.de_max, g.de_max);

    // airspeed -> throttle
    u.dt = std::clamp(trim.dt + g.k_V * (Va_c - Va), 0.0, 1.0);
}

void LongitudinalAutopilot::reset() {
    h_int = 0.0;
    theta_c = theta_trim;
}

} // namespace BSIM
