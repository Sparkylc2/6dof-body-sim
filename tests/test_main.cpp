#include "autopilot.hpp"
#include "body.hpp"
#include "body_system.hpp"
#include "control.hpp"
#include "force_models/fixed_wing_model.hpp"
#include "frames.hpp"
#include "sensors.hpp"
#include "world.hpp"
#include <cmath>
#include <cstdio>
#include <functional>
#include <memory>
#include <numbers>
#include <random>

// moving forward i'm having claude generate the tests, in the interest of time
// however, i do go through the tests, and some of them are my own suggestions
static int g_checks = 0;
static int g_failures = 0;

// these i had claude generate for me
#define CHECK(cond)                                                                                \
    do {                                                                                           \
        ++g_checks;                                                                                \
        if (!(cond)) {                                                                             \
            ++g_failures;                                                                          \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);                            \
        }                                                                                          \
    } while (0)

#define CHECK_NEAR(a, b, tol)                                                                      \
    do {                                                                                           \
        ++g_checks;                                                                                \
        const double a_ = (a), b_ = (b);                                                           \
        if (!(std::abs(a_ - b_) <= (tol))) {                                                       \
            ++g_failures;                                                                          \
            std::printf("FAIL %s:%d: %s = %.17g, expected %.17g\n", __FILE__, __LINE__, #a, a_,    \
                        b_);                                                                       \
        }                                                                                          \
    } while (0)

#define CHECK_VEC_NEAR(a, b, tol)                                                                  \
    do {                                                                                           \
        ++g_checks;                                                                                \
        const Eigen::Vector3d a_ = (a), b_ = (b);                                                  \
        if (!((a_ - b_).norm() <= (tol))) {                                                        \
            ++g_failures;                                                                          \
            std::printf("FAIL %s:%d: %s = [%.12g %.12g %.12g], expected [%.12g %.12g %.12g]\n",    \
                        __FILE__, __LINE__, #a, a_.x(), a_.y(), a_.z(), b_.x(), b_.y(), b_.z());   \
        }                                                                                          \
    } while (0)

constexpr double pi = std::numbers::pi;

// the full inertia matrix from the packed [Ix, Iy, Iz, Ixz]
// with Ixy = Iyz = 0 from symmetry about the x-z plane
static Eigen::Matrix3d inertia_matrix(const Eigen::Vector4d &I) {
    Eigen::Matrix3d J;
    J << I(0), 0, -I(3), 0, I(1), 0, -I(3), 0, I(2);
    return J;
}

// a body with the zagi's mass properties whose forces and moments are set by each test
// (no gravity unless a test adds it)
struct TestBody : uav::Body {
    std::function<Eigen::Vector3d(const uav::Body &, const uav::BodyState &)> get_forces;
    std::function<Eigen::Vector3d(const uav::Body &, const uav::BodyState &)> get_moments;

    TestBody() {
        m = 1.56;
        I << 0.1147, 0.0576, 0.1712, 0.0015;
    }

    uav::Wrench wrench(const uav::BodyState &x, const uav::AirState &) const override {
        uav::Wrench w;
        if (get_forces)
            w.F = get_forces(*this, x);
        if (get_moments)
            w.M = get_moments(*this, x);
        return w;
    }
};

// a system holding one zagi at rest, level, at the origin
// tests read and write the body through *sys.bodies[0]
static uav::BodySystem make_zagi_system() {
    uav::BodySystem sys;
    sys.bodies.push_back(std::make_unique<TestBody>());
    return sys;
}

// the test body inside a system made by make_zagi_system
static TestBody &test_body(uav::BodySystem &sys) { return static_cast<TestBody &>(*sys.bodies[0]); }

// ---- frames (euler) ----
static void test_rotation_is_orthonormal() {
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> ang(-pi, pi);
    double worst_orth = 0.0, worst_det = 0.0;
    for (int i = 0; i < 1000; ++i) {
        const Eigen::Matrix3d R = uav::rotation_euler(ang(rng), ang(rng) / 2.0, ang(rng));
        worst_orth = std::max(worst_orth, (R * R.transpose() - Eigen::Matrix3d::Identity()).norm());
        worst_det = std::max(worst_det, std::abs(R.determinant() - 1.0));
    }
    CHECK_NEAR(worst_orth, 0.0, 1e-12);
    CHECK_NEAR(worst_det, 0.0, 1e-12);
}

// rotation_euler() is R_v^b (NED -> body) so the body axes in NED are the columns of the transpose
static void test_rotation_sign_conventions() {
    const Eigen::Vector3d nose = uav::rotation_euler(0, 0, pi / 2).transpose().col(0);
    CHECK_VEC_NEAR(nose, Eigen::Vector3d(0, 1, 0), 1e-12); // psi = 90 deg means nose points east

    const Eigen::Vector3d nose_pitched = uav::rotation_euler(0, pi / 6, 0).transpose().col(0);
    CHECK(nose_pitched.z() < 0.0); // theta > 0 means nose up (NED d is down)

    const Eigen::Vector3d right_wing = uav::rotation_euler(pi / 6, 0, 0).transpose().col(1);
    CHECK(right_wing.z() > 0.0); // phi > 0 meaning the right wing down
}

// at phi = theta = 0 the Euler rates equal the body rates
static void test_euler_rates_identity_when_level() {
    const Eigen::Matrix3d E = uav::body_rates_to_euler_rates(0.0, 0.0);
    CHECK_NEAR((E - Eigen::Matrix3d::Identity()).norm(), 0.0, 1e-15);
}

// ---- frames (quaternion) ----
// rotation(e) is R_b^v, so for the same attitude it has to be the transpose of eqn 2.5
static void test_quat_rotation_matches_euler() {
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> ang(-pi, pi);
    double worst = 0.0;
    for (int i = 0; i < 1000; ++i) {
        const double phi = ang(rng), theta = ang(rng) / 2.0, psi = ang(rng);
        const Eigen::Matrix3d R_q = uav::rotation(uav::euler_to_quat(phi, theta, psi));
        const Eigen::Matrix3d R_e = uav::rotation_euler(phi, theta, psi).transpose();
        worst = std::max(worst, (R_q - R_e).norm());
    }
    CHECK_NEAR(worst, 0.0, 1e-12);
}

// euler -> quat -> euler gives back the same angles (away from theta = +-90 deg)
static void test_quat_euler_round_trip() {
    std::mt19937 rng(11);
    std::uniform_real_distribution<double> ang(-pi, pi);
    double worst = 0.0;
    for (int i = 0; i < 1000; ++i) {
        const Eigen::Vector3d euler(ang(rng), 0.99 * ang(rng) / 2.0, ang(rng));
        const Eigen::Vector3d back =
            uav::quat_to_euler(uav::euler_to_quat(euler.x(), euler.y(), euler.z()));
        worst = std::max(worst, (back - euler).norm());
    }
    CHECK_NEAR(worst, 0.0, 1e-9);
}

// quat_rate (B.5) has to agree with eqn 3.16: moving the euler angles at the 3.16 rates should
// move the quaternion at quat_rate (central difference)
static void test_quat_rate_matches_euler_rates() {
    std::mt19937 rng(3);
    std::uniform_real_distribution<double> ang(-pi, pi);
    std::uniform_real_distribution<double> rate(-2.0, 2.0);
    double worst = 0.0;
    for (int i = 0; i < 1000; ++i) {
        const Eigen::Vector3d euler(ang(rng), 0.8 * ang(rng) / 2.0, ang(rng));
        const Eigen::Vector3d w(rate(rng), rate(rng), rate(rng));
        const Eigen::Vector3d euler_dot = uav::body_rates_to_euler_rates(euler.x(), euler.y()) * w;

        const double h = 1e-6;
        const Eigen::Vector3d ep = euler + h * euler_dot, em = euler - h * euler_dot;
        const Eigen::Vector4d fd = (uav::euler_to_quat(ep.x(), ep.y(), ep.z()) -
                                    uav::euler_to_quat(em.x(), em.y(), em.z())) /
                                   (2.0 * h);
        const Eigen::Vector4d e = uav::euler_to_quat(euler.x(), euler.y(), euler.z());
        worst = std::max(worst, (fd - uav::quat_rate(e, w)).norm());
    }
    CHECK_NEAR(worst, 0.0, 1e-8);
}

// ---- dynamics ----
// if u have zero body rates, v_dot = f / m
static void test_constant_force_gives_linear_velocity() {
    auto sys = make_zagi_system();
    const double a = 2.0;
    test_body(sys).get_forces = [a](const uav::Body &body, const uav::BodyState &) {
        return Eigen::Vector3d(body.m * a, 0, 0);
    };

    const double dt = 1e-3;
    for (int i = 0; i < 1000; ++i)
        sys.step(dt);

    const uav::Body &b = *sys.bodies[0];
    CHECK_VEC_NEAR(b.state.vel(), Eigen::Vector3d(a * 1.0, 0, 0), 1e-12); // u = a t after 1 s
    CHECK_NEAR(b.state.pos().x(), 0.5 * a, 1e-9); // x = a t^2 / 2 (rk4 is exact for this)
    CHECK_NEAR(b.state.pos().y(), 0.0, 1e-15);
    CHECK_NEAR(b.state.pos().z(), 0.0, 1e-15);
}

// position rates are R_b^v [u v w]^T, and its yawed 90 deg and flying forward
// along the body x axis, the craft must move east
static void test_yawed_forward_flight_moves_east() {
    auto sys = make_zagi_system();
    uav::Body &b = *sys.bodies[0];
    b.state.quat() = uav::euler_to_quat(0, 0, pi / 2);
    b.state.vel() = Eigen::Vector3d(10.0, 0, 0);

    const double dt = 1e-3;
    for (int i = 0; i < 1000; ++i)
        sys.step(dt);

    CHECK_VEC_NEAR(b.state.pos(), Eigen::Vector3d(0, 10.0, 0), 1e-9); // 10 m east after 1 s
}

// if you start from rest the gyro terms should vanish, so w_dot = I^-1*tau.
// this checks gamma, gamma_3, gamma_4, gamma_8 (eqn. 3.13) and the 1/Iy term against
// I.inverse()
static void test_angular_accel_from_rest_matches_Jinv() {
    auto sys = make_zagi_system();
    const uav::Body &b = *sys.bodies[0];
    const Eigen::Vector3d moment(0.3, -0.2, 0.1);

    const uav::StateVector x_dot = sys.derivatives(b, b.state, Eigen::Vector3d::Zero(), moment);

    const Eigen::Vector3d expected = inertia_matrix(b.I).inverse() * moment;
    CHECK_VEC_NEAR(x_dot.segment<3>(uav::IDX::RATES), expected, 1e-12);
}

// eqn 3.10 is I * w_dot = m - w x (I*w).
// with no applied moment w_dot should equal I^-1 (-w x I*w).
// this checks gamma_1, gamma_2, gamma_5, gamma_6, gamma_7 (eqn. 3.13)
static void test_gyroscopic_terms_match_vector_form() {
    auto sys = make_zagi_system();
    uav::Body &b = *sys.bodies[0];
    const Eigen::Vector3d w0(0.5, -0.3, 0.8);
    b.state.rates() = w0;

    const uav::StateVector x_dot =
        sys.derivatives(b, b.state, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());

    const Eigen::Matrix3d J = inertia_matrix(b.I);
    const Eigen::Vector3d expected = J.inverse() * (-w0.cross(J * w0));
    CHECK_VEC_NEAR(x_dot.segment<3>(uav::IDX::RATES), expected, 1e-12);
}

// v_dot = -w x v (eqn 3.15) with no forces
static void test_coriolis_term_matches_vector_form() {
    auto sys = make_zagi_system();
    uav::Body &b = *sys.bodies[0];
    const Eigen::Vector3d w0(0.5, -0.3, 0.8), v0(12.0, 1.0, -2.0);
    b.state.rates() = w0;
    b.state.vel() = v0;

    const uav::StateVector x_dot =
        sys.derivatives(b, b.state, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());

    CHECK_VEC_NEAR(x_dot.segment<3>(uav::IDX::VEL), -w0.cross(v0), 1e-12);
}

// with Ixz 0 a pure roll moment shouldn't give any yaw rates
static void test_roll_yaw_coupling() {
    for (const double Ixz : {0.0, 0.0015}) {
        auto sys = make_zagi_system();
        TestBody &b = test_body(sys);
        b.I.w() = Ixz;
        b.get_moments = [](const uav::Body &, const uav::BodyState &) {
            return Eigen::Vector3d(0.1, 0, 0);
        };

        for (int i = 0; i < 100; ++i)
            sys.step(1e-3);

        if (Ixz == 0.0)
            CHECK(b.state.rates().z() == 0.0);
        else
            CHECK(std::abs(b.state.rates().z()) > 1e-6);
    }
}

// a torque free body keeps its rotational energy 1/2 w^T I w and the size of its angular momentum
// |I w|, so this checks the gyro terms and the rk4 together over a long run
static void test_torque_free_conserves_energy_and_momentum() {
    auto sys = make_zagi_system();
    uav::Body &b = *sys.bodies[0];
    b.state.rates() = Eigen::Vector3d(0.5, -0.3, 0.8);

    const Eigen::Matrix3d J = inertia_matrix(b.I);
    auto energy = [&J](const Eigen::Vector3d &w) { return 0.5 * w.dot(J * w); };
    const double E0 = energy(b.state.rates()), H0 = (J * b.state.rates()).norm();

    for (int i = 0; i < 10000; ++i)
        sys.step(1e-3);

    CHECK_NEAR(energy(b.state.rates()) / E0, 1.0, 1e-9);
    CHECK_NEAR((J * b.state.rates()).norm() / H0, 1.0, 1e-9);
}

// the attitude quaternion stays unit length while tumbling
static void test_quat_stays_unit() {
    auto sys = make_zagi_system();
    uav::Body &b = *sys.bodies[0];
    b.state.rates() = Eigen::Vector3d(1.0, -2.0, 0.5);

    for (int i = 0; i < 5000; ++i)
        sys.step(1e-2);

    CHECK_NEAR(b.state.quat().norm(), 1.0, 1e-12);
}

// the whole point of the quaternions, pitching straight up through theta = 90 deg (where 3.16
// blows up) has to work. with only q nonzero the gyro terms are all zero so q stays constant,
// and after pi/2 seconds at 1 rad/s the nose points straight up (-d in NED)
static void test_pitch_through_vertical() {
    auto sys = make_zagi_system();
    uav::Body &b = *sys.bodies[0];
    b.state.rates() = Eigen::Vector3d(0, 1.0, 0);

    const int n = 1000;
    const double dt = (pi / 2.0) / n;
    for (int i = 0; i < n; ++i)
        sys.step(dt);

    const Eigen::Vector3d nose = uav::rotation(b.state.quat()).col(0);
    CHECK_VEC_NEAR(nose, Eigen::Vector3d(0, 0, -1), 1e-9);

    // and keep going past vertical, it should end up upside down facing south after pi seconds
    for (int i = 0; i < n; ++i)
        sys.step(dt);

    const Eigen::Vector3d nose_over = uav::rotation(b.state.quat()).col(0);
    CHECK_VEC_NEAR(nose_over, Eigen::Vector3d(-1, 0, 0), 1e-9);
}

// ---- bodies ----
// with no keyboard input the drone only feels gravity, straight down the body z axis when level
// and along -x when pitched straight up
static void test_drone_gravity_only() {
    const uav::Drone d;

    const uav::Wrench level = d.wrench(uav::BodyState(), {});
    CHECK_VEC_NEAR(level.F, Eigen::Vector3d(0, 0, d.m * 9.81), 1e-12);
    CHECK_VEC_NEAR(level.M, Eigen::Vector3d::Zero(), 1e-15);

    uav::BodyState pitched;
    pitched.quat() = uav::euler_to_quat(0, pi / 2, 0);
    CHECK_VEC_NEAR(d.wrench(pitched, {}).F, Eigen::Vector3d(-d.m * 9.81, 0, 0), 1e-12);
}

// ---- fixed wing force models ----
// air data in still air with the body level, so the ground velocity is the air velocity
static uav::AirData still_air_data(const Eigen::Vector3d &vel_b) {
    return uav::air_data(vel_b, Eigen::Vector4d(1, 0, 0, 0), uav::AirState());
}

// forward along x gives no alpha or beta, tilting the velocity down body z gives alpha, and at rest
// everything is 0 (not nan)
static void test_air_data() {
    const uav::AirData fwd = still_air_data(Eigen::Vector3d(10, 0, 0));
    CHECK_NEAR(fwd.Va, 10.0, 1e-12);
    CHECK_NEAR(fwd.alpha, 0.0, 1e-15);
    CHECK_NEAR(fwd.beta, 0.0, 1e-15);

    const uav::AirData pitched = still_air_data(Eigen::Vector3d(10, 0, 10 * std::tan(0.1)));
    CHECK_NEAR(pitched.alpha, 0.1, 1e-12);

    const uav::AirData slipping = still_air_data(Eigen::Vector3d(10, 10 * std::tan(0.2), 0));
    CHECK_NEAR(slipping.beta, 0.2, 1e-12);

    const uav::AirData rest = still_air_data(Eigen::Vector3d::Zero());
    CHECK(rest.Va == 0.0 && rest.alpha == 0.0 && rest.beta == 0.0);
}

// well below stall the blended lift is the linear lift, past stall it drops off
static void test_lift_stall() {
    const uav::FixedWingParams p = uav::generate_zagi_flying_wing().params;

    for (const double alpha : {-0.1, 0.0, 0.05, 0.15})
        CHECK_NEAR(uav::lift_coefficient(p, alpha), p.CL0 + p.CL_alpha * alpha, 1e-6);

    CHECK(uav::lift_coefficient(p, 0.7) < uav::lift_coefficient(p, 0.3));
}

// well below stall the blended drag is the polar, well past stall it's the flat plate drag
static void test_drag_stall() {
    const uav::FixedWingParams p = uav::generate_zagi_flying_wing().params;
    const double AR = p.b * p.b / p.S;

    for (const double alpha : {-0.1, 0.0, 0.05, 0.15}) {
        const double CL = p.CL0 + p.CL_alpha * alpha;
        CHECK_NEAR(uav::drag_coefficient(p, alpha), p.CD_p + CL * CL / (pi * p.e * AR), 1e-6);
    }

    for (const double alpha : {-1.2, 1.2, pi / 2}) {
        const double s = std::sin(alpha);
        CHECK_NEAR(uav::drag_coefficient(p, alpha), 2.0 * s * s * std::abs(s), 1e-6);
    }
}

// longitudinal force at q = de = 0 has to be lift and drag rotated by alpha
static void test_longitudinal_is_rotated_lift_drag() {
    const uav::FixedWingParams p = uav::generate_zagi_flying_wing().params;
    const double alpha = 0.12;
    const uav::AirData air{15.0, alpha, 0.0};

    const double qbar_S = 0.5 * air.rho * air.Va * air.Va * p.S;
    const double L = qbar_S * uav::lift_coefficient(p, alpha);
    const double D = qbar_S * uav::drag_coefficient(p, alpha);

    const uav::Wrench w = uav::longitudinal_aero(p, air, 0.0, 0.0);
    CHECK_NEAR(w.F.x(), -D * std::cos(alpha) + L * std::sin(alpha), 1e-12);
    CHECK_NEAR(w.F.z(), -D * std::sin(alpha) - L * std::cos(alpha), 1e-12);
    CHECK(w.F.z() < 0.0); // lift points up (-d)
}

// cm_alpha < 0 so more alpha should pitch the nose back down (static stability)
static void test_pitch_stability() {
    const uav::FixedWingParams p = uav::generate_zagi_flying_wing().params;
    const double m0 = uav::longitudinal_aero(p, {15.0, 0.0, 0.0}, 0.0, 0.0).M.y();
    const double m1 = uav::longitudinal_aero(p, {15.0, 0.1, 0.0}, 0.0, 0.0).M.y();
    CHECK(m1 < m0);
}

// with no sideslip, roll or yaw rate, aileron or rudder there's nothing lateral
static void test_lateral_symmetry() {
    const uav::FixedWingParams p = uav::generate_zagi_flying_wing().params;
    const uav::Wrench w = uav::lateral_aero(p, {15.0, 0.1, 0.0}, Eigen::Vector3d(0, 0.3, 0), 0, 0);
    CHECK_VEC_NEAR(w.F, Eigen::Vector3d::Zero(), 1e-15);
    CHECK_VEC_NEAR(w.M, Eigen::Vector3d::Zero(), 1e-15);
}

// static thrust is 1/2 rho S_prop C_prop k_motor^2 dt^2 and it hits zero when Va = k_motor dt
static void test_propulsion() {
    const uav::FixedWingParams p = uav::generate_zagi_flying_wing().params;
    const uav::AirData rest{0.0, 0.0, 0.0};
    CHECK_NEAR(uav::propulsion(p, rest, 1.0).F.x(),
               0.5 * rest.rho * p.S_prop * p.C_prop * p.k_motor * p.k_motor, 1e-12);
    CHECK_NEAR(uav::propulsion(p, {0.5 * p.k_motor, 0.0, 0.0}, 0.5).F.x(), 0.0, 1e-12);
}

// at rest (even spinning) the zagi only feels gravity, and nothing is nan
static void test_zagi_at_rest_is_gravity_only() {
    uav::FixedWing zagi = uav::generate_zagi_flying_wing();
    zagi.state.rates() = Eigen::Vector3d(1.0, -0.5, 0.3);
    zagi.controls = {0.1, -0.1, 0.0, 0.0};

    const uav::Wrench w = zagi.wrench(zagi.state, {});
    CHECK_VEC_NEAR(w.F, Eigen::Vector3d(0, 0, zagi.m * 9.81), 1e-12);
    CHECK_VEC_NEAR(w.M, Eigen::Vector3d::Zero(), 1e-12);
}

// flying it through the rk4 for a few seconds shouldn't blow up
static void test_zagi_flies_without_nan() {
    uav::BodySystem sys;
    auto zagi = std::make_unique<uav::FixedWing>(uav::generate_zagi_flying_wing());
    zagi->state.vel() = Eigen::Vector3d(15.0, 0, 0);
    zagi->controls.dt = 0.8;
    uav::FixedWing &b = *zagi;
    sys.bodies.push_back(std::move(zagi));

    for (int i = 0; i < 3000; ++i)
        sys.step(1e-3);

    CHECK(b.state.state.allFinite());
}

// ---- wind ----
// the wind is taken off the ground velocity, steady wind rotated into the body and gusts as they
// are
static void test_air_data_with_wind() {
    // flying north at 10 m/s ground speed into a 5 m/s headwind (wind blowing south) is 15 m/s
    // airspeed
    uav::AirState head;
    head.wind_i = Eigen::Vector3d(-5, 0, 0);
    const Eigen::Vector4d level(1, 0, 0, 0);
    CHECK_NEAR(uav::air_data(Eigen::Vector3d(10, 0, 0), level, head).Va, 15.0, 1e-12);

    // yawed 90 deg (nose east) the same wind comes from the left (north), so relative to the air
    // the body is slipping left (v_r < 0) and beta < 0
    const Eigen::Vector4d east = uav::euler_to_quat(0, 0, pi / 2);
    const uav::AirData cross = uav::air_data(Eigen::Vector3d(10, 0, 0), east, head);
    CHECK_NEAR(cross.Va, std::sqrt(125.0), 1e-12);
    CHECK_NEAR(cross.beta, -std::atan2(5.0, 10.0), 1e-12);

    // an upward gust in body z (w_g < 0) gives the wing more alpha
    uav::AirState gust;
    gust.gust_b = Eigen::Vector3d(0, 0, -1.0);
    CHECK_NEAR(uav::air_data(Eigen::Vector3d(10, 0, 0), level, gust).alpha, std::atan2(1.0, 10.0),
               1e-12);

    // density is passed straight through
    uav::AirState thin;
    thin.rho = 1.0;
    CHECK_NEAR(uav::air_data(Eigen::Vector3d(10, 0, 0), level, thin).rho, 1.0, 0.0);
}

// the standard deviation of each gust channel over a long run
static Eigen::Vector3d gust_std(double Va, double dt, double T) {
    uav::DrydenGust g;
    std::mt19937 rng(1);
    const int n = static_cast<int>(T / dt);
    Eigen::Vector3d sum = Eigen::Vector3d::Zero(), sum_sq = Eigen::Vector3d::Zero();
    for (int i = 0; i < n; ++i) {
        const Eigen::Vector3d x = g.step(Va, dt, rng);
        sum += x;
        sum_sq += x.cwiseProduct(x);
    }
    const Eigen::Vector3d mean = sum / n;
    return (sum_sq / n - mean.cwiseProduct(mean)).cwiseSqrt();
}

// the dryden filters are normalised so each gust comes out with rms sigma, and the sqrt(dt) noise
// scaling means that doesn't change with the step size
// (10000 s is about 1000 correlation times for u, so the estimate is good to a few percent)
static void test_dryden_gust_rms() {
    const uav::DrydenParams p;
    const Eigen::Vector3d sigma(p.sigma_u, p.sigma_v, p.sigma_w);
    for (const double dt : {0.02, 0.01})
        CHECK_VEC_NEAR(gust_std(20.0, dt, 10000.0).cwiseQuotient(sigma), Eigen::Vector3d::Ones(),
                       0.08);
}

// with turbulence off and no wind the world has to fly the exact same trajectory as the bodies
// stepping on their own
static void test_world_calm_matches_body_system() {
    auto make_zagi = []() {
        auto z = std::make_unique<uav::FixedWing>(uav::generate_zagi_flying_wing());
        z->state.vel() = Eigen::Vector3d(15.0, 0, 0);
        z->controls = {-0.1, 0.02, 0.0, 0.8};
        return z;
    };

    uav::BodySystem sys;
    sys.bodies.push_back(make_zagi());
    uav::World world;
    world.atmosphere.turbulence = false;
    world.bodies.bodies.push_back(make_zagi());

    for (int i = 0; i < 2000; ++i) {
        sys.step(1e-3);
        world.step(1e-3);
    }

    CHECK((sys.bodies[0]->state.state - world.bodies.bodies[0]->state.state).norm() == 0.0);
}

// with turbulence on the world still flies without nans, and actually gets bumped around
static void test_world_turbulent_flight() {
    uav::World world;
    auto zagi = std::make_unique<uav::FixedWing>(uav::generate_zagi_flying_wing());
    zagi->state.vel() = Eigen::Vector3d(15.0, 0, 0);
    zagi->controls.dt = 0.8;
    uav::FixedWing &b = *zagi;
    world.bodies.bodies.push_back(std::move(zagi));

    for (int i = 0; i < 3000; ++i)
        world.step(1e-3);

    CHECK(b.state.state.allFinite());
    CHECK(b.air.gust_b.norm() > 0.0);
}

// ---- sensors ----
// all the noise and biases off so the readings should be exactly the truth
static uav::SensorParams ideal_sensors() {
    uav::SensorParams p;
    p.sigma_gyro = p.sigma_accel = p.sigma_abs_pres = p.sigma_diff_pres = p.sigma_mag = 0.0;
    p.beta_abs_pres = p.beta_diff_pres = p.beta_mag = 0.0;
    p.sigma_gps = Eigen::Vector3d::Zero();
    p.sigma_gps_V = 0.0;
    return p;
}

// mean and std of a stream of samples
struct Stats {
    double sum = 0.0, sum_sq = 0.0;
    int n = 0;
    void add(double x) {
        sum += x;
        sum_sq += x * x;
        ++n;
    }
    double mean() const { return sum / n; }
    double std() const { return std::sqrt(sum_sq / n - mean() * mean()); }
};

// 7.3 (forces) needs to give the same specific force as 7.1 (state derivatives)
// this js uses a zagi with some bank, pitch, and rotation with some sideslip
static void test_specific_force_matches_state_derivatives() {
    uav::FixedWing zagi = uav::generate_zagi_flying_wing();
    zagi.state.quat() = uav::euler_to_quat(0.4, 0.15, -0.7);
    zagi.state.vel() = Eigen::Vector3d(17.0, 1.5, 1.2);
    zagi.state.rates() = Eigen::Vector3d(0.3, -0.2, 0.15);
    zagi.controls = {-0.1, 0.05, 0.0, 0.8};

    const uav::BodyState &x = zagi.state;
    const uav::Wrench w = zagi.wrench(x, {});
    uav::BodySystem sys;
    const uav::StateVector x_dot = sys.derivatives(zagi, x, w.F, w.M);

    const double u = x.vel().x(), v = x.vel().y(), wv = x.vel().z();
    const double p = x.rates().x(), q = x.rates().y(), r = x.rates().z();
    const Eigen::Vector3d euler = uav::quat_to_euler(x.quat());
    const double phi = euler.x(), theta = euler.y(), g = uav::GRAVITY;

    // eqn 7.1 without the noise
    const Eigen::Vector3d a(
        x_dot(uav::VEL) + q * wv - r * v + g * std::sin(theta),
        x_dot(uav::VEL + 1) + r * u - p * wv - g * std::cos(theta) * std::sin(phi),
        x_dot(uav::VEL + 2) + p * v - q * u - g * std::cos(theta) * std::cos(phi));

    CHECK_VEC_NEAR(uav::specific_force(w, zagi.m, x.quat()), a, 1e-10);
}

// in free fall the accelerometer shld read nothing, but if it's held up against gravity it should
// read -g on the axis pointing down
static void test_specific_force_free_fall_and_hover() {
    uav::Drone d;
    const Eigen::Vector4d level(1, 0, 0, 0);
    CHECK_VEC_NEAR(uav::specific_force(uav::gravity(d.m, level), d.m, level),
                   Eigen::Vector3d::Zero(), 1e-12);

    // zero net force (eg sitting on a table) is the same reading for any attitude, -R^T [0 0 g]
    const Eigen::Vector4d tilted = uav::euler_to_quat(0.3, -0.2, 1.0);
    const Eigen::Vector3d expect = -uav::rotation(tilted).transpose() * Eigen::Vector3d(0, 0, 9.81);
    CHECK_VEC_NEAR(uav::specific_force(uav::Wrench(), d.m, tilted), expect, 1e-12);
    CHECK_VEC_NEAR(uav::specific_force(uav::Wrench(), d.m, level), Eigen::Vector3d(0, 0, -9.81),
                   1e-12);
}

// simple calcs for what the pressure should be given airspeed, rho and alt
static void test_pressure_models() {
    CHECK_NEAR(uav::abs_pressure(1.2682, 100.0), 1.2682 * 9.81 * 100.0, 1e-9);
    CHECK_NEAR(uav::abs_pressure(1.2682, 0.0), 0.0, 0.0);
    CHECK_NEAR(uav::diff_pressure(1.2682, 20.0), 0.5 * 1.2682 * 400.0, 1e-9);
}

// rotating provo's magnetic field (7.4) into the body and tilt compensating it has to give back
// psi - declination for any roll and pitch
static void test_magnetic_heading_tilt_compensation() {
    const Eigen::Vector3d m_i(21053.0, 4520.0, 47689.0); // nT, north east down
    const double declination = std::atan2(m_i.y(), m_i.x());
    CHECK_NEAR(declination, 12.12 * uav::DEG, 0.01 * uav::DEG);

    std::mt19937 rng(3);
    std::uniform_real_distribution<double> angle(-1.2, 1.2);
    // check over a couple different tilt angles
    for (int i = 0; i < 50; ++i) {
        const double phi = angle(rng), theta = angle(rng), psi = 2.5 * angle(rng);
        const Eigen::Vector3d m_b =
            uav::rotation(uav::euler_to_quat(phi, theta, psi)).transpose() * m_i;
        const double psi_m = uav::magnetic_heading(m_b, phi, theta);
        CHECK_NEAR(uav::wrap_angle(psi_m + declination - psi), 0.0, 1e-12);
    }
}

// the gauss-markov process should settle to sigma / sqrt(1 - e^(-2*k*Ts)), and with table 7.2 that
// should be the 4.7 m and 9.2 m bias sizes from 7.5.2 (the book quotes 9.4 m for altitude)
static void test_gauss_markov() {
    CHECK_NEAR(uav::gauss_markov_std(1.0 / 1100.0, 1.0, 0.21), 4.7, 0.3);
    CHECK_NEAR(uav::gauss_markov_std(1.0 / 1100.0, 1.0, 0.40), 9.4, 0.05);

    // a shorter time constant  (rather than 1100s xD) so a long run covers lots of correlation
    // times
    const double k = 0.1, Ts = 1.0, sigma = 0.4;
    std::mt19937 rng(5);
    Stats s;
    double nu = 0.0;
    for (int i = 0; i < 200000; ++i) {
        nu = uav::gauss_markov(nu, k, Ts, sigma, rng);
        s.add(nu);
    }
    CHECK_NEAR(s.std() / uav::gauss_markov_std(k, Ts, sigma), 1.0, 0.03);
    CHECK_NEAR(s.mean(), 0.0, 0.1);
}

static void test_ground_speed_course() {
    const Eigen::Vector2d east = uav::ground_speed_course(Eigen::Vector3d(0, 15, -2));
    CHECK_NEAR(east.x(), 15.0, 1e-12); // climbing doesn't count
    CHECK_NEAR(east.y(), pi / 2, 1e-12);
    const Eigen::Vector2d sw = uav::ground_speed_course(Eigen::Vector3d(-3, -4, 0));
    CHECK_NEAR(sw.x(), 5.0, 1e-12);
    CHECK_NEAR(sw.y(), std::atan2(-4.0, -3.0), 1e-12);
}

// a zagi with an empty sensor object in the world
static uav::FixedWing &add_zagi(uav::World &world, const uav::SensorParams &params) {
    auto zagi = std::make_unique<uav::FixedWing>(uav::generate_zagi_flying_wing());
    zagi->state.pos() = Eigen::Vector3d(10, -5, -30);
    zagi->state.quat() = uav::euler_to_quat(0.2, 0.05, pi / 2); // nose east
    zagi->state.vel() = Eigen::Vector3d(15.0, 0.5, 1.0);
    zagi->controls = {-0.1, 0.0, 0.0, 0.8};
    zagi->sensors.params = params;
    uav::FixedWing &b = *zagi;
    world.bodies.bodies.push_back(std::move(zagi));
    return b;
}

// with the noise off every reading should be the model applied to the true state
static void test_ideal_sensors_read_truth() {
    uav::World world;
    world.atmosphere.turbulence = false;
    world.atmosphere.wind_i = Eigen::Vector3d(-3, 2, 0);
    uav::FixedWing &b = add_zagi(world, ideal_sensors());

    world.step(0.01);

    const uav::BodyState &x = b.state;
    const uav::SensorReadings &y = b.sensors.y;
    CHECK_VEC_NEAR(y.gyro, x.rates(), 1e-15);
    CHECK_VEC_NEAR(y.accel, uav::specific_force(b.wrench(x, b.air), b.m, x.quat()), 1e-12);
    CHECK_NEAR(y.abs_pres, b.air.rho * 9.81 * -x.pos().z(), 1e-9);
    const double Va = uav::air_data(x.vel(), x.quat(), b.air).Va;
    CHECK_NEAR(y.diff_pres, 0.5 * b.air.rho * Va * Va, 1e-9);
    CHECK_NEAR(y.compass, uav::quat_to_euler(x.quat()).z(), 1e-12);
    CHECK_VEC_NEAR(y.gps_pos, Eigen::Vector3d(x.pos().x(), x.pos().y(), -x.pos().z()), 1e-12);

    const Eigen::Vector3d vel_i = uav::rotation(x.quat()) * x.vel();
    CHECK_NEAR(y.gps_Vg, vel_i.head<2>().norm(), 1e-12);
    CHECK_NEAR(y.gps_course, std::atan2(vel_i.y(), vel_i.x()), 1e-12);
    CHECK(std::abs(y.gps_course - pi / 2) < 0.2); // heading east so roughly east

    // the wind changes the airspeed (and so diff_pres) but not what gps sees, which is ground speed
    CHECK(std::abs(Va - x.vel().norm()) > 1.0);
}

// each sensor group samples at its own rate
// for analog 2s at 100Hz is 200 analog samples
// the compass is at 8Hz, it samples on the first step then every 0.125s
// the gps is at 1Hz, it samples on the first step then every 1s
static void test_sensor_sample_rates() {
    uav::World world;
    world.atmosphere.turbulence = false;
    uav::FixedWing &b = add_zagi(world, uav::SensorParams());

    int n_analog = 0, n_compass = 0, n_gps = 0;
    uav::SensorReadings last = b.sensors.y;
    for (int i = 0; i < 200; ++i) {
        world.step(0.01);
        const uav::SensorReadings &y = b.sensors.y;
        n_analog += y.t_analog != last.t_analog;
        n_compass += y.t_compass != last.t_compass;
        n_gps += y.t_gps != last.t_gps;

        // an unsampled reading shouldn't change between samples
        if (y.t_gps == last.t_gps)
            CHECK(y.gps_pos == last.gps_pos);
        last = y;
    }
    CHECK(n_analog == 200);
    CHECK(n_compass == 17);
    CHECK(n_gps == 3);
}

// a drone hovering still at 20 m
// the readings should scatter around with truth + bias from the appendix H standard deviations
static void test_sensor_noise_statistics() {
    uav::World world;
    world.atmosphere.turbulence = false;
    auto drone = std::make_unique<uav::Drone>();
    drone->state.pos() = Eigen::Vector3d(0, 0, -20);
    drone->state.quat() = uav::euler_to_quat(0, 0, 0.5);
    drone->input.F = Eigen::Vector3d(0, 0, -drone->m * 9.81); // holds it up still
    uav::Drone &d = *drone;
    world.bodies.bodies.push_back(std::move(drone));

    const uav::SensorParams p;
    Stats gyro, accel_z, abs_pres, diff_pres, compass;
    double t_compass = -1.0;
    for (int i = 0; i < 20000; ++i) {
        world.step(0.01);
        const uav::SensorReadings &y = d.sensors.y;
        gyro.add(y.gyro.y());
        accel_z.add(y.accel.z());
        abs_pres.add(y.abs_pres);
        diff_pres.add(y.diff_pres);
        if (y.t_compass != t_compass) {
            compass.add(y.compass);
            t_compass = y.t_compass;
        }
    }

    // makes sure it really did stay still
    CHECK_VEC_NEAR(d.state.pos(), Eigen::Vector3d(0, 0, -20), 1e-9);

    CHECK_NEAR(gyro.mean(), 0.0, 0.05 * p.sigma_gyro);
    CHECK_NEAR(gyro.std() / p.sigma_gyro, 1.0, 0.05);
    CHECK_NEAR(accel_z.mean(), -9.81, 0.05 * p.sigma_accel);
    CHECK_NEAR(accel_z.std() / p.sigma_accel, 1.0, 0.05);
    CHECK_NEAR(abs_pres.mean(), d.air.rho * 9.81 * 20.0 + p.beta_abs_pres, 0.05 * p.sigma_abs_pres);
    CHECK_NEAR(abs_pres.std() / p.sigma_abs_pres, 1.0, 0.05);
    CHECK_NEAR(diff_pres.mean(), p.beta_diff_pres, 0.05 * p.sigma_diff_pres); // still air, Va = 0
    CHECK_NEAR(diff_pres.std() / p.sigma_diff_pres, 1.0, 0.05);
    CHECK_NEAR(compass.mean(), 0.5 + p.beta_mag, 0.1 * p.sigma_mag);
    CHECK_NEAR(compass.std() / p.sigma_mag, 1.0, 0.1);

    // 200 s of gps is a fraction of the 1100 s time constant, so the error is mostly the bias it
    // started with, so this js checks its a sensible size
    CHECK((d.sensors.y.gps_pos - Eigen::Vector3d(0, 0, 20)).norm() < 50.0);
}

// ---- control ----
static void test_saturate_and_rate_limit() {
    CHECK_NEAR(uav::saturate(2.0, 0.0, 1.0), 1.0, 0.0);
    CHECK_NEAR(uav::saturate(-0.5, 0.0, 1.0), 0.0, 0.0);
    CHECK_NEAR(uav::saturate(0.3, 0.0, 1.0), 0.3, 0.0);

    // 2 per second over 0.1 s can only move 0.2 either way
    CHECK_NEAR(uav::rate_limit(1.0, 0.0, 2.0, 0.1), 0.2, 1e-15);
    CHECK_NEAR(uav::rate_limit(-1.0, 0.0, 2.0, 0.1), -0.2, 1e-15);
    CHECK_NEAR(uav::rate_limit(0.05, 0.0, 2.0, 0.1), 0.05, 1e-15);
}

// a 100 Hz clock should sample 100 times a second whether the sim steps at Ts or faster
// (over [0, 1) since the first sample is at t = 0)
static void test_sample_clock() {
    for (const double dt : {0.01, 0.001, 0.004}) {
        uav::SampleClock clock{0.01};
        int n = 0;
        const int steps = static_cast<int>(std::round(1.0 / dt));
        for (int i = 0; i < steps; ++i)
            n += clock.ready(i * dt);
        CHECK(n == 100);
    }
}

// the trapezoid rule is exact for a ramp, and the derivative should converge to the slope
static void test_integrator_and_differentiator() {
    const double dt = 0.01;
    uav::Integrator integ;
    uav::Differentiator diff;
    for (int i = 0; i <= 100; ++i) {
        integ.step(3.0 * i * dt, dt);
        diff.step(3.0 * i * dt, dt);
    }
    CHECK_NEAR(integ.y, 1.5, 1e-12); // integral of 3t from 0 to 1
    CHECK_NEAR(diff.y, 3.0, 1e-6);   // 1 s is 20 time constants of the filter

    // an angle going round through +-pi at 1 rad/s should still differentiate to 1
    uav::Differentiator ang;
    ang.wrap = true;
    for (int i = 0; i <= 200; ++i)
        ang.step(uav::wrap_angle(2.5 + i * dt), dt);
    CHECK_NEAR(ang.y, 1.0, 1e-6);
}

// a low pass with a step gets to 1 - 1/e after one time constant, and should be 1 at steady state
static void test_low_pass() {
    uav::LowPass lp;
    lp.tau = 0.5;
    const double dt = 0.001;
    lp.step(0.0, dt); // primes at 0 so the step comes after
    for (int i = 0; i < 500; ++i)
        lp.step(1.0, dt);
    CHECK_NEAR(lp.y, 1.0 - std::exp(-1.0), 1e-3);
    for (int i = 0; i < 10000; ++i)
        lp.step(1.0, dt);
    CHECK_NEAR(lp.y, 1.0, 1e-9);
}

// commanding 179 deg while at -179 deg should be a 2 deg error not 358
static void test_pid_wraps_angles() {
    uav::PID pid;
    pid.params = {.kp = 1.0, .wrap = true};
    CHECK_NEAR(pid.update(179.0 * uav::DEG, -179.0 * uav::DEG, 0.01), -2.0 * uav::DEG, 1e-12);
}

// a big step into a saturated pi on an pure integrator, without anti-windup the integrator
// should build up the entire time
static double windup_overshoot(uav::AntiWindup aw) {
    uav::PID pid;
    pid.params = {.kp = 1.0, .ki = 1.0, .lo = -1.0, .hi = 1.0, .anti_windup = aw, .Tt = 0.2};
    double x = 0.0, x_max = 0.0;
    const double dt = 0.01;
    // 120 s since without anti-windup it takes a long time to unwind and settle
    for (int i = 0; i < 12000; ++i) {
        x += dt * pid.update(10.0, x, dt);
        x_max = std::max(x_max, x);
    }
    CHECK_NEAR(x, 10.0, 1e-3); // all of them get there eventually (we hope)
    return x_max - 10.0;
}

static void test_pid_anti_windup() {
    const double none = windup_overshoot(uav::AntiWindup::None);
    const double clamp = windup_overshoot(uav::AntiWindup::Clamp);
    CHECK(none > 5.0);
    CHECK(clamp < 0.5);
}

// flying north at 15 m/s ground speed through a crosswind from the east, the nose points into the
// wind but the course should still be north
static void test_true_control_state() {
    uav::BodyState x;
    x.pos() = Eigen::Vector3d(5, -3, -40);
    const double crab = std::atan2(4.0, 15.0); // nose turned east into the wind
    x.quat() = uav::euler_to_quat(0, 0, crab);
    x.vel() = uav::rotation(x.quat()).transpose() * Eigen::Vector3d(15, 0, 0);
    uav::AirState air;
    air.wind_i = Eigen::Vector3d(0, -4, 0); // blowing west

    const uav::ControlState s = uav::true_control_state(x, air);
    CHECK_NEAR(s.h, 40.0, 1e-12);
    CHECK_NEAR(s.chi, 0.0, 1e-12);
    CHECK_NEAR(s.Vg, 15.0, 1e-12);
    CHECK_NEAR(s.psi, crab, 1e-12);
    CHECK_NEAR(s.Va, std::sqrt(15.0 * 15.0 + 16.0), 1e-12);
}

// the autopilot runs on its clock and holds control outputs in between
struct CountingAutopilot : uav::Autopilot {
    int n = 0;
    uav::FixedWingControls update(const uav::AutopilotCommands &,
                                  const uav::ControlState &) override {
        ++n;
        return {.dt = 0.01 * n};
    }
    void reset() override { n = 0; }
};

static void test_autopilot_holds_between_samples() {
    CountingAutopilot ap;
    ap.clock.Ts = 0.05;
    for (int i = 0; i < 100; ++i)
        ap.step(i * 0.01, {}, {});
    CHECK(ap.n == 20);
    CHECK_NEAR(ap.controls.dt, 0.2, 1e-12);
}

int main() {
    test_rotation_is_orthonormal();
    test_rotation_sign_conventions();
    test_euler_rates_identity_when_level();
    test_quat_rotation_matches_euler();
    test_quat_euler_round_trip();
    test_quat_rate_matches_euler_rates();
    test_constant_force_gives_linear_velocity();
    test_yawed_forward_flight_moves_east();
    test_angular_accel_from_rest_matches_Jinv();
    test_gyroscopic_terms_match_vector_form();
    test_coriolis_term_matches_vector_form();
    test_roll_yaw_coupling();
    test_torque_free_conserves_energy_and_momentum();
    test_quat_stays_unit();
    test_pitch_through_vertical();
    test_drone_gravity_only();
    test_air_data();
    test_lift_stall();
    test_drag_stall();
    test_longitudinal_is_rotated_lift_drag();
    test_pitch_stability();
    test_lateral_symmetry();
    test_propulsion();
    test_zagi_at_rest_is_gravity_only();
    test_zagi_flies_without_nan();
    test_air_data_with_wind();
    test_dryden_gust_rms();
    test_world_calm_matches_body_system();
    test_world_turbulent_flight();
    test_specific_force_matches_state_derivatives();
    test_specific_force_free_fall_and_hover();
    test_pressure_models();
    test_magnetic_heading_tilt_compensation();
    test_gauss_markov();
    test_ground_speed_course();
    test_ideal_sensors_read_truth();
    test_sensor_sample_rates();
    test_sensor_noise_statistics();
    test_saturate_and_rate_limit();
    test_sample_clock();
    test_integrator_and_differentiator();
    test_low_pass();
    test_pid_wraps_angles();
    test_pid_anti_windup();
    test_true_control_state();
    test_autopilot_holds_between_samples();

    std::printf("%d/%d passed\n (pray its all right xD)", g_checks - g_failures, g_checks);
    return g_failures == 0 ? 0 : 1;
}
