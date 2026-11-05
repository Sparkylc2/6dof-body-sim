#include "sensors.hpp"
#include "body.hpp"
#include "force_models/fixed_wing_model.hpp"
#include "frames.hpp"
#include <algorithm>
#include <cmath>

namespace BSIM {

Vector3d specific_force(const Wrench &w, double m, const Vector4d &e) {
    // the mass doesn't feel gravity relative to the case (it pulls on both the same),
    // so the accelerometer only reads (F - F_gravity)/m

    // gravity() already gives force in body frame
    return (w.F - gravity(m, e).F) / m;
}

double abs_pressure(double rho, double h_agl) {
    // P_ground - P = rho*g*h_agl
    return rho * GRAVITY * h_agl;
}

double diff_pressure(double rho, double Va) {
    // P_t - P_s = 1/2*rho*Va^2
    return 0.5 * rho * Va * Va;
}

double magnetic_heading(const Vector3d &m_b, double phi, double theta) {
    const double c_ph = std::cos(phi), s_ph = std::sin(phi);
    const double c_th = std::cos(theta), s_th = std::sin(theta);

    // R_b^v1 = R_v2^v1(theta) R_b^v2(phi), takes the body measurement back to the vehicle-1 frame
    // which is only yawed, so the x-y part of it is horizontal
    // from eqn 7.12
    const Eigen::Matrix3d R{
        {c_th, s_th * s_ph, s_th * c_ph},
        {0.0, c_ph, -s_ph},
        {-s_th, c_th * s_ph, c_th * c_ph},
    };
    const Vector3d m_v1 = R * m_b;

    // the field sits at -psi_m in the v1 frame, eqn 7.13
    return -std::atan2(m_v1.y(), m_v1.x());
}

double gauss_markov(double nu, double k, double Ts, double sigma, std::mt19937 &rng) {
    // standard gauss markov model
    std::normal_distribution<double> noise(0.0, 1.0);
    return std::exp(-k * Ts) * nu + sigma * noise(rng);
}

double gauss_markov_std(double k, double Ts, double sigma) {
    // var(.)= a^2*var(.)+ sigma^2 at steady state with a = e^(-k*Ts)
    return sigma / std::sqrt(1.0 - std::exp(-2.0 * k * Ts));
}

Vector2d ground_speed_course(const Vector3d &vel_i) {
    const double Vn = vel_i.x(), Ve = vel_i.y();
    // atan2 instead of tan^-1 so the course comes out in the right quadrant
    return {std::sqrt(Vn * Vn + Ve * Ve), std::atan2(Ve, Vn)};
}

double wrap_angle(double a) { return std::remainder(a, 2.0 * std::numbers::pi); }

void Sensors::update(const Body &b, double t, std::mt19937 &rng) {
    std::normal_distribution<double> noise(0.0, 1.0);
    auto noise3 = [&]() { return Vector3d(noise(rng), noise(rng), noise(rng)); };

    const SensorParams &p = params;
    const AircraftState &x = b.state;

    // a sensor is ready to be measured once t gets to next sample time,
    // if the sim step is longer than the sensor period the missed samples are just dropped
    // (the 1e-9 is so floating point build up in next doesn't make a sample that's exactly due look
    // a hair early, when the sim steps at exactly Ts that sample would otherwise be lost)
    auto ready = [t](double &next, double Ts) {
        if (t + 1e-9 < next)
            return false;

        next += Ts;
        if (next <= t)
            next = t + Ts;

        return true;
    };

    // analog sensors, they get sampled every Ts
    if (ready(next_analog, p.Ts)) {
        // the rate gyros, eqn 7.5
        const Vector3d gyro = x.rates() + p.bias_gyro + p.sigma_gyro * noise3();
        // the gyros have a built in limit
        y.gyro = gyro.cwiseMax(-p.range_gyro).cwiseMin(p.range_gyro);

        // accelerometers, the forces are already computed by the sim so the state
        // derivatives for vel aren't needed (unlike in 7.1)
        const Vector3d f = specific_force(b.wrench(x, b.air), b.m, x.quat());
        const Vector3d accel = f + p.bias_accel + p.sigma_accel * noise3();
        y.accel = accel.cwiseMax(-p.range_accel).cwiseMin(p.range_accel);

        // absolute pressure, eqn 7.9, the ground is at p_d = 0 so h_agl = -p_d
        y.abs_pres =
            abs_pressure(b.air.rho, -x.pos().z()) + p.beta_abs_pres + p.sigma_abs_pres * noise(rng);

        // diff. pressures, eqn 7.10, Va is relative to the air so wind and gusts neeed to be
        // accounted for
        const AirData air = air_data(x.vel(), x.quat(), b.air);
        y.diff_pres =
            diff_pressure(air.rho, air.Va) + p.beta_diff_pres + p.sigma_diff_pres * noise(rng);

        y.t_analog = t;
    }

    // digital compass, gets sampled every Ts_compass
    if (ready(next_compass, p.Ts_compass)) {
        // eqn 7.14, the true heading with a declination bias and magnetometer noise
        const double psi = quat_to_euler(x.quat()).z();
        y.compass = wrap_angle(psi + p.beta_mag + p.sigma_mag * noise(rng));
        y.t_compass = t;
    }

    // gps, every Ts_gps
    if (ready(next_gps, p.Ts_gps)) {
        for (int i = 0; i < 3; ++i) {
            if (y.t_gps < 0.0) {
                // the time constant is like 1100s so i'm having it start somewhere on its steady
                // state distribution
                nu_gps(i) = gauss_markov_std(p.k_gps, p.Ts_gps, p.sigma_gps(i)) * noise(rng);
            } else {
                nu_gps(i) = gauss_markov(nu_gps(i), p.k_gps, p.Ts_gps, p.sigma_gps(i), rng);
            }
        }

        // position, eqn 7.18 - 7.20, h is up so it's -p_d
        y.gps_pos = Vector3d(x.pos().x(), x.pos().y(), -x.pos().z()) + nu_gps;

        // ground speed and course, eqn 7.25 and 7.26
        // the sim already has the actual ground velocities so this uses those rather than the
        // no sideslip/no climb formula the textbook uses
        const Vector2d Vg_chi = ground_speed_course(rotation(x.quat()) * x.vel());

        // eqn 7.24, the course reading gets worse with 1/Vg, so it's clamped so it's not infinite
        // when the body has no velocity
        const double sigma_chi = p.sigma_gps_V / std::max(Vg_chi.x(), 1.0);
        y.gps_Vg = Vg_chi.x() + p.sigma_gps_V * noise(rng);
        y.gps_course = wrap_angle(Vg_chi.y() + sigma_chi * noise(rng));

        y.t_gps = t;
    }
}

} // namespace BSIM
