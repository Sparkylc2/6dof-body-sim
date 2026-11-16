#pragma once
#include <limits>

namespace BSIM {

constexpr double INF = std::numeric_limits<double>::infinity();

// --- basic blocks ---
// clamps u into [lo, hi]
double saturate(double u, double lo, double hi);

// limits how fast a signal changes between samples (a servo can only move so fast)
// the rate is per second
double rate_limit(double u, double u_prev, double rate, double dt);

// keeps track of when a discrete thing (controller, filter etc) is ready to be run in the sim
// if dt is longer than Ts the missed samples are just dropped
struct SampleClock {
    double Ts = 0.01;
    double next = 0.0;

    // returns true (and updates next) if a sample is ready at time t
    bool ready(double t);
};

// --- discrete filters ---
// all use the trapezoid rule which keeps it stable for any dt and matches
// the continuous version well as long as dt is small compared to the time constants
//
// the first sample after a reset just readies the previous input so there's no jump from 0
// trapezoidal integrator, y[n] = y[n-1] + dt/2 (u[n] + u[n-1])
struct Integrator {
    double y = 0.0;
    double u_prev = 0.0;
    bool ready = false;

    double step(double u, double dt);
    void reset(double y0 = 0.0);
};

// band limited differentiator s / (tau s + 1),
// a pure s isn't causal and amplifies noise so tau fixes that above 1/tau rad/s
// (wrap is for when the input is an angle, so going through +-pi isn't seen as a huge jump)
// [queue flight dynamics autopilot lab flashbacks]
struct Differentiator {
    double tau = 0.05; // todo: play around w this, this is just a guess taken from the lab
    bool wrap = false;

    double y = 0.0;
    double u_prev = 0.0;
    bool ready = false;

    double step(double u, double dt);
    void reset();
};

// first order low pass 1 / (tau s + 1), with unity dc gain
struct LowPass {
    double tau = 0.1;

    double y = 0.0;
    double u_prev = 0.0;
    bool ready = false;

    double step(double u, double dt);
    void reset(double y0 = 0.0);
};

// --- pid ---
// sets what the derivative term acts on. iirc you use measurement to prevent a step in setpoint
// from causing a derivative kick but investigate later.  error is just kd*de/dt
enum class DerivativeOn { Error };

// how the integrator is stopped from winding up while the output is saturated
// clamp stops integrating when saturated
enum class AntiWindup { None, Clamp };

// parameters for a PID controller
struct PIDParams {
    double kp = 0.0, ki = 0.0, kd = 0.0;
    double lo = -INF, hi = INF; // output limits

    double tau = 0.05; // derivative filter time constant
    DerivativeOn derivative = DerivativeOn::Error;

    AntiWindup anti_windup = AntiWindup::Clamp;
    double Tt = 0.1; // tracking time constant for back calculation

    bool wrap = false; // the error is an angle so it gets wrapped into [-pi, pi]
};

// a discrete pid
//   u = ff + kp*e + integral(ki*e) + kd (d/dt), saturated into [lo, hi]
// e = r - y, ff is any feedforward the controller receives (eg a trim value)
//
// the integrator holds the integral of ki*e (not ki times the integral of e)
// so changing ki mid flight shouldn't cause it to jump
struct PID {
    PIDParams params;

    Integrator integral;  // integral of e
    Differentiator deriv; // derivative of e (for now)

    double u = 0.0; // last output

    // differentiates internally (filtered)
    double update(double r, double y, double dt, double ff = 0.0);

    // clears the integrator and filters
    void reset();

private:
    double error(double r, double y) const;
    // shared by both updates once the derivative term is known
    double output(double e, double d_term, double dt, double ff);
};

} // namespace BSIM
