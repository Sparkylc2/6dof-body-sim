#include "lateral_autopilot.hpp"
#include "force_models.hpp" // GRAVITY
#include <cmath>

namespace BSIM {

namespace {
double sign(double x) { return (x > 0.0) - (x < 0.0); }
} // namespace

LateralTransferFunctions lateral_transfer_functions(const FixedWingParams &p, double m,
                                                    const Vector4d &J, double rho, double Va) {
    // roll moment coefficients mixed with the yaw ones through the inertia (the Gamma terms)
    const double Jx = J(0), Jz = J(2), Jxz = J(3);
    const double g1 = Jx * Jz - Jxz * Jxz;
    const double g3 = Jz / g1;
    const double g4 = Jxz / g1;
    const double Cp_p = g3 * p.Cl_p + g4 * p.Cn_p;
    const double Cp_da = g3 * p.Cl_da + g4 * p.Cn_da;

    LateralTransferFunctions tf;
    // eqn 5.23 and 5.24
    tf.a_phi1 = -0.5 * rho * Va * Va * p.S * p.b * Cp_p * p.b / (2.0 * Va);
    tf.a_phi2 = 0.5 * rho * Va * Va * p.S * p.b * Cp_da;
    // sideslip model, section 5.4
    tf.a_beta1 = -rho * Va * p.S / (2.0 * m) * p.CY_beta;
    tf.a_beta2 = rho * Va * p.S / (2.0 * m) * p.CY_dr;
    return tf;
}

RollGains roll_gains(const LateralTransferFunctions &tf, const LateralDesignParams &d) {
    RollGains g;
    g.kp = d.da_max / d.e_phi_max * sign(tf.a_phi2);                // eqn 6.7
    g.wn = std::sqrt(std::abs(tf.a_phi2) * d.da_max / d.e_phi_max); // eqn 6.8
    g.kd = (2.0 * d.zeta_phi * g.wn - tf.a_phi1) / tf.a_phi2;       // eqn 6.9
    g.ki = d.ki_phi;
    return g;
}

CourseGains course_gains(const RollGains &roll, const LateralDesignParams &d) {
    CourseGains g;
    g.wn = roll.wn / d.W_chi;
    g.kp = 2.0 * d.zeta_chi * g.wn * d.Vg / GRAVITY; // eqn 6.12
    g.ki = g.wn * g.wn * d.Vg / GRAVITY;             // eqn 6.13
    return g;
}

SideslipGains sideslip_gains(const LateralTransferFunctions &tf, const LateralDesignParams &d) {
    SideslipGains g;
    g.kp = d.dr_max / d.e_beta_max * sign(tf.a_beta2);             // eqn 6.16
    g.wn = (tf.a_beta1 + tf.a_beta2 * g.kp) / (2.0 * d.zeta_beta); // from eqn 6.15
    g.ki = g.wn * g.wn / tf.a_beta2;                               // eqn 6.17 (via 6.14)
    return g;
}

// --- hold loops ---
void RollHold::set_gains(const RollGains &g, double da_max) {
    pid.params.kp = g.kp;
    pid.params.ki = g.ki;
    pid.params.kd = 0.0; // damping comes from p directly, see update
    pid.params.lo = -da_max;
    pid.params.hi = da_max;
    kd = g.kd;
}

double RollHold::update(double phi_c, double phi, double p, double dt) {
    // -kd p goes in as the feedforward so it's inside the saturation and anti windup stuff
    return pid.update(phi_c, phi, dt, -kd * p);
}

void CourseHold::set_gains(const CourseGains &g, double phi_max) {
    pid.params.kp = g.kp;
    pid.params.ki = g.ki;
    pid.params.kd = 0.0;
    pid.params.lo = -phi_max;
    pid.params.hi = phi_max;
    pid.params.wrap = true;
}

double CourseHold::update(double chi_c, double chi, double dt) {
    return pid.update(chi_c, chi, dt);
}

void SideslipHold::set_gains(const SideslipGains &g, double dr_max) {
    pid.params.kp = g.kp;
    pid.params.ki = g.ki;
    pid.params.kd = 0.0;
    pid.params.lo = -dr_max;
    pid.params.hi = dr_max;
}

double SideslipHold::update(double beta, double dt) {
    // beta_c = 0, so the error becomes -beta
    return pid.update(0.0, beta, dt);
}

// --- lateral autopilot ---
void LateralAutopilot::design(const FixedWingParams &p, double m, const Vector4d &J, double rho,
                              double Va, const LateralDesignParams &d) {
    tf = lateral_transfer_functions(p, m, J, rho, Va);

    roll_g = roll_gains(tf, d);
    course_g = course_gains(roll_g, d);
    roll.set_gains(roll_g, d.da_max);
    course.set_gains(course_g, d.phi_max);

    has_rudder = tf.a_beta2 != 0.0;
    if (has_rudder) {
        sideslip_g = sideslip_gains(tf, d);
        sideslip.set_gains(sideslip_g, d.dr_max);
    }
}

void LateralAutopilot::step(double t, double chi_c, double phi, double chi, double beta, double p,
                            FixedWingControls &u) {
    if (!clock.ready(t)) {
        return;
    }
    const double dt = clock.Ts;
    phi_c = course.update(chi_c, chi, dt);
    u.da = roll.update(phi_c, phi, p, dt);
    if (has_rudder) {
        u.dr = sideslip.update(beta, dt);
    }
}

void LateralAutopilot::reset() {
    roll.reset();
    course.reset();
    sideslip.reset();
    phi_c = 0.0;
}

} // namespace BSIM
