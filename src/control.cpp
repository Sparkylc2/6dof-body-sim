#include "control.hpp"
#include "sensors.hpp" // wrap_angle
#include <algorithm>

namespace uav {

// --- basic blocks ---
double saturate(double u, double lo, double hi) { return std::clamp(u, lo, hi); }

double rate_limit(double u, double u_prev, double rate, double dt) {
    const double max_step = rate * dt;
    return u_prev + std::clamp(u - u_prev, -max_step, max_step);
}

bool SampleClock::ready(double t) {
    next += Ts;
    // if it's fallen behind by more than a sample it just skips
    if (next <= t) {
        next = t + Ts;
    }
    return true;
}

// --- discrete filters ---
double Integrator::step(double u, double dt) {
    if (!ready) {
        u_prev = u;
        ready = true;
    }
    // area of the trapezoid under u over the step
    y += dt / 2.0 * (u + u_prev);
    u_prev = u;
    return y;
}

void Integrator::reset(double y0) {
    y = y0;
    u_prev = 0.0;
    ready = false;
}

double Differentiator::step(double u, double dt) {
    if (!ready) {
        u_prev = u;
        ready = true;
    }
    // s becomes (2/dt) (z - 1)/(z + 1) (s in s / (tau s + 1)) and rearranging
    // y[n] = (2 tau - dt)/(2 tau + dt) y[n-1] + 2/(2 tau + dt) (u[n] - u[n-1])
    const double du = wrap ? wrap_angle(u - u_prev) : u - u_prev;
    y = (2.0 * tau - dt) / (2.0 * tau + dt) * y + 2.0 / (2.0 * tau + dt) * du;
    u_prev = u;
    return y;
}

void Differentiator::reset() {
    y = 0.0;
    u_prev = 0.0;
    ready = false;
}

double LowPass::step(double u, double dt) {
    if (!ready) {
        u_prev = u;
        ready = true;
    }
    // same trapezoidal substitution into 1 / (tau s + 1)
    // y[n] = (2 tau - dt)/(2 tau + dt) y[n-1] + dt/(2 tau + dt) (u[n] + u[n-1])
    y = (2.0 * tau - dt) / (2.0 * tau + dt) * y + dt / (2.0 * tau + dt) * (u + u_prev);
    u_prev = u;
    return y;
}

void LowPass::reset(double y0) {
    y = y0;
    u_prev = 0.0;
    ready = false;
}

// --- pid ---
double PID::error(double r, double y) const {
    const double e = r - y;
    return params.wrap ? wrap_angle(e) : e;
}

double PID::update(double r, double y, double dt, double ff) {
    const PIDParams &p = params;
    const double e = error(r, y);

    // the filter settings are in params so they can be changed in one place
    deriv.tau = p.tau;
    deriv.wrap = p.wrap;

    const double d_term = p.kd * deriv.step(e, dt);
    return output(e, d_term, dt, ff);
}

double PID::output(double e, double d_term, double dt, double ff) {
    const PIDParams &p = params;

    const double I_prev = integral.y;
    integral.step(p.ki * e, dt);

    // basic pid formula
    const double u_unsat = ff + p.kp * e + integral.y + d_term;
    u = saturate(u_unsat, p.lo, p.hi);

    switch (p.anti_windup) {
    case AntiWindup::None:
        break;

    case AntiWindup::Clamp:
        // only undoes this step's integration if it's gone further past the limit
        // if it hasn't it can start unwinding immediately
        if ((u_unsat > p.hi && p.ki * e > 0.0) || (u_unsat < p.lo && p.ki * e < 0.0)) {
            integral.y = I_prev;
            u = saturate(ff + p.kp * e + integral.y + d_term, p.lo, p.hi);
        }
        break;
    }

    return u;
}

void PID::reset() {
    integral.reset();
    deriv.reset();
    u = 0.0;
}

} // namespace uav
