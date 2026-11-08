#pragma once
#include "body_state.hpp"
#include "control.hpp"
#include "force_models/fixed_wing_model.hpp"

namespace BSIM {

// --- lateral autopilot (chapter 6.3) ---
// the roll attitude hold (inner loop), course hold (outer loop, commands roll) and sideslip hold
// gains come from the closed form design equations, so they js need the transfer function
// coefficients and the tuning parameters
struct LateralTransferFunctions {
    double a_phi1 = 0.0, a_phi2 = 0.0;
    double a_beta1 = 0.0, a_beta2 = 0.0;
};

// J is the inertia [Jx, Jy, Jz, Jxz], rho and Va are density and airspeed
// from eqns 5.23, 5.24 (roll) and the sideslip model in section 5.4
LateralTransferFunctions lateral_transfer_functions(const FixedWingParams &p, double m,
                                                    const Vector4d &J, double rho, double Va);

// --- tuning parameters ---
struct LateralDesignParams {
    // roll
    double da_max = 0.0;    // aileron limit
    double e_phi_max = 0.0; // roll error that saturates the ailerons
    double zeta_phi = 0.0;  // roll damping ratio
    double ki_phi = 0.0;    // roll integrator gain (from the root locus, not closed form)

    // course
    double W_chi = 0.0;    // bandwidth sep from the roll loop, w_n_phi = W_chi w_n_chi
    double zeta_chi = 0.0; // course damping ratio
    double Vg = 0.0;       // nominal ground speed
    double phi_max = 0.0;  // roll angle limit

    // sideslip
    double dr_max = 0.0;     // rudder limit
    double e_beta_max = 0.0; // sideslip error that saturates the rudder
    double zeta_beta = 0.0;  // sideslip damping ratio
};

// --- gains ---
struct RollGains {
    double kp = 0.0, ki = 0.0, kd = 0.0;
    double wn = 0.0; // roll loop nat freq
};

struct CourseGains {
    double kp = 0.0, ki = 0.0;
    double wn = 0.0; // course loop nat freq
};

struct SideslipGains {
    double kp = 0.0, ki = 0.0;
    double wn = 0.0; // sideslip nat freq
};

// eqns 6.7 (kp), 6.8 (wn) and 6.9 (kd), ki gets passed through
RollGains roll_gains(const LateralTransferFunctions &tf, const LateralDesignParams &d);

// w_n_chi = w_n_phi / W_chi, then eqns 6.12 (kp) and 6.13 (ki)
CourseGains course_gains(const RollGains &roll, const LateralDesignParams &d);

// eqns 6.16 (kp) and 6.17 (ki), wn from 6.15
SideslipGains sideslip_gains(const LateralTransferFunctions &tf, const LateralDesignParams &d);

// --- hold loops ---
// each one is a PID with the output limited to the actuator/roll limit

// da = kp (phi_c - phi) + ki/s (phi_c - phi) - kd p
// in this one the derivative acts on the measured roll rate instead of error
struct RollHold {
    PID pid;
    double kd = 0.0;

    void set_gains(const RollGains &g, double da_max);
    double update(double phi_c, double phi, double p, double dt);
    void reset() { pid.reset(); }
};

// phi_c = kp (chi_c - chi) + ki/s (chi_c - chi)
// and the error is wrapped into [-pi, pi]
struct CourseHold {
    PID pid;

    void set_gains(const CourseGains &g, double phi_max);
    double update(double chi_c, double chi, double dt);
    void reset() { pid.reset(); }
};

// dr = -kp beta - ki/s beta
// (also clear to see it gives no input when beta is 0)
struct SideslipHold {
    PID pid;

    void set_gains(const SideslipGains &g, double dr_max);
    double update(double beta, double dt);
    void reset() { pid.reset(); }
};

// --- lateral autopilot ---
// course hold to roll hold to ailerons, and sideslip hold to rudder
// only sets da and dr
struct LateralAutopilot {
    SampleClock clock; // clock.Ts is the dt for the loops

    RollHold roll;
    CourseHold course;
    SideslipHold sideslip;
    bool has_rudder = false;

    LateralTransferFunctions tf;
    RollGains roll_g;
    CourseGains course_g;
    SideslipGains sideslip_g;

    double phi_c = 0.0; // last roll cmd

    // works out the transfer functions at (rho, Va) and sets the gains based on that
    void design(const FixedWingParams &p, double m, const Vector4d &J, double rho, double Va,
                const LateralDesignParams &d);

    // if a sample is ready at t, updates da (and dr if there's a rudder) in u
    // chi_c is the commanded course, the rest is what the autopilot receives
    // (either true state or estimate)
    void step(double t, double chi_c, double phi, double chi, double beta, double p,
              FixedWingControls &u);

    void reset();
};

} // namespace BSIM
