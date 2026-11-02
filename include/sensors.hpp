#pragma once
#include <Eigen/Dense>
#include <numbers>
#include <random>

#include "force_models.hpp"

namespace uav {
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;

// the sensors need the whole body to sample, but the body holds the sensor
struct Body;

constexpr double DEG = std::numbers::pi / 180.0;

// --- sensor parameters ---
// noise, bias and sample rates for all the sensors in chapter 7,
// the defaults are from the example parts from appendix H
// everything is SI (rad, m/s^2, Pa, m)
struct SensorParams {
    // gyros, accels and the pressure sensors are analog, they get sampled through the
    // control loop (50 to 100 Hz), the compass and gps are digital (and slower)
    double Ts = 0.01;          // control loop
    double Ts_compass = 0.125; // HMR3300 runs at 8 Hz (H.4)
    double Ts_gps = 1.0;       // H.5

    // rate gyros, ADXRS450 (H.1), sigma = noise density * sqrt(bandwidth)
    double sigma_gyro = 0.13 * DEG;
    Vector3d bias_gyro = Vector3d::Zero(); // 7.5 assumes the bias is calibrated out
    double range_gyro = 350.0 * DEG;

    // accelerometers, ADXL325 (H.2), sigma = (250 mu * g)[noise density]/sqrt(Hz) * sqrt(bandwidth)
    double sigma_accel = 0.0025 * GRAVITY;
    Vector3d bias_accel = Vector3d::Zero(); // same logic as the gyros
    double range_accel = 6.0 * GRAVITY;

    // pressure sensors (H.3), the betas are the temperature drift that's left after calibration
    double beta_abs_pres = 125.0, sigma_abs_pres = 10.0; // MP3H6115A, absolute (altitude)
    double beta_diff_pres = 20.0, sigma_diff_pres = 2.0; // MPXV5004G, differential (airspeed)

    // digital compass, HMR3300 (H.4), beta is mostly from not knowing the declination exactly
    double beta_mag = 1.0 * DEG, sigma_mag = 0.3 * DEG;

    // gps position, gauss-markov parameters from table 7.2
    double k_gps = 1.0 / 1100.0;          // 1 / time constant of the error
    Vector3d sigma_gps{0.21, 0.21, 0.40}; // std of the driving noise for [n, e, h]
    // gps velocity, 7.5.3 gives 0.01 to 0.05 m/s per axis so i'm just using the worst end
    double sigma_gps_V = 0.05;
};

// autopilot gets fed this, each reading stays until the sensor samples again
struct SensorReadings {
    Vector3d gyro = Vector3d::Zero();    // [p, q, r] (7.5)
    Vector3d accel = Vector3d::Zero();   // specific force [ax, ay, az] in the body frame (7.3)
    double abs_pres = 0.0;               // P_ground - P (7.9)
    double diff_pres = 0.0;              // P_t - P_s (7.10)
    double compass = 0.0;                // heading psi (7.14)
    Vector3d gps_pos = Vector3d::Zero(); // [n, e, h] (7.18 - 7.20)
    double gps_Vg = 0.0;                 // ground speed (7.25)
    double gps_course = 0.0;             // course chi (7.26)

    // when each group of sensors was last sampled (should be useful for the state estimation)
    // (analog is the gyros, accels and pressure sensors), and -1 means never
    double t_analog = -1.0, t_compass = -1.0, t_gps = -1.0;
};

// a body's full sensor suite
struct Sensors {
    SensorParams params;
    SensorReadings y;

    Vector3d nu_gps = Vector3d::Zero(); // gauss-markov gps position error [nu_n, nu_e, nu_h]

    // when each group is next ready to be sampled
    double next_analog = 0.0, next_compass = 0.0, next_gps = 0.0;

    // samples whatever sensors are ready at time t from the body at that time
    void update(const Body &b, double t, std::mt19937 &rng);
};

// --- noise free sensor models ---
// the accelerometer feels everything but gravity, 7.1 (in the body frame)
Vector3d specific_force(const Wrench &w, double m, const Vector4d &e);

// pressure drop from the ground up to h_agl assuming the density doesn't change
// from eqn 7.7 (and 7.9) (todo: add isa here instead, should be quick)
double abs_pressure(double rho, double h_agl);

// dynamic pressure (ie what a pitot tube would measure)
// from bernoulli in eqn 7.10
double diff_pressure(double rho, double Va);

// magnetic heading from the body magnetometers (accounting for roll and pitch, so js levelling out)
// from eqn 7.12 and 7.13 (todo: doesn't include declination, could add that later [eqn 7.11])
double magnetic_heading(const Vector3d &m_b, double phi, double theta);

// one sample of the gauss-markov process nu[n+1] = e^(-k*Ts) nu[n] + eta, eta ~ N(0, sigma^2)
// from eqn 7.17, (used for gps error)
double gauss_markov(double nu, double k, double Ts, double sigma, std::mt19937 &rng);

// the steady state std the gauss-markov process arrives at, sigma / sqrt(1 - e^(-2*k*Ts))
double gauss_markov_std(double k, double Ts, double sigma);

// ground speed and course [Vg, chi] from the inertial ground velocity
// from eqn 7.21 and 7.22
Vector2d ground_speed_course(const Vector3d &vel_i);

// wraps an angle into [-pi, pi]
double wrap_angle(double a);

} // namespace uav
