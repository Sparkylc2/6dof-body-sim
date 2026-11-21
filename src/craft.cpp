#include "craft.hpp"
#include <cmath>
#include <numbers>

namespace BSIM {

// --- drone ---
Drone::Drone() {
    m = 1.56;

    const double Jx = 0.1147, Jy = 0.0576, Jz = 0.1712, Jxz = 0.0015;

    // following the symmetry assumption that airplanes use
    J << Jx, Jy, Jz, Jxz;
}

Wrench Drone::wrench(const AircraftState &x, const AirState &) const {
    Wrench w = gravity(m, x.quat());
    w += input;
    return w;
}

// --- fixed wing ---
FixedWing::FixedWing(double m, const Vector4d &J, const FixedWingParams &params) : params(params) {
    this->m = m;
    this->J = J;
}

Wrench FixedWing::wrench(const AircraftState &x, const AirState &air) const {
    if (aileron.ideal()) {
        return wrench(x, air, controls);
    }
    FixedWingControls u = controls;
    u.da = aileron.pos;
    return wrench(x, air, u);
}

void FixedWing::step_actuators(double dt) {
    if (!aileron.ideal()) {
        aileron.step(controls.da, dt);
    }
}

Wrench FixedWing::wrench(const AircraftState &x, const AirState &air,
                         const FixedWingControls &u) const {
    return fixed_wing_wrench(m, params, x, air, u);
}

// factories for the fixed wing
FixedWing generate_zagi_flying_wing() {
    const double m = 1.56;
    const double Jx = 0.1147, Jy = 0.0576, Jz = 0.1712, Jxz = 0.0015;

    // the zagi is a flying wing with elevons and no rudder, so all the dr coefficients stay 0
    const FixedWingParams p{
        .S = 0.2589,
        .b = 1.4224,
        .c = 0.3302,
        .e = 0.9,

        .CL0 = 0.09167,
        .CL_alpha = 3.5016,
        .CL_q = 2.8932,
        .CL_de = 0.2724,
        .CD0 = 0.01631,
        .CD_alpha = 0.2108,
        .CD_q = 0.0,
        .CD_de = 0.3045,
        .CD_p = 0.0254,
        .Cm0 = -0.02338,
        .Cm_alpha = -0.5675,
        .Cm_q = -1.3990,
        .Cm_de = -0.3254,

        .M = 50.0,
        .alpha0 = 0.4712,

        .CY0 = 0.0,
        .CY_beta = -0.07359,
        .CY_p = 0.0,
        .CY_r = 0.0,
        .CY_da = 0.0,
        .Cl0 = 0.0,
        .Cl_beta = -0.02854,
        .Cl_p = -0.3209,
        .Cl_r = 0.03066,
        .Cl_da = 0.1682,
        .Cn0 = 0.0,
        .Cn_beta = -0.00040,
        .Cn_p = -0.01297,
        .Cn_r = -0.00434,
        .Cn_da = -0.00328,

        .S_prop = 0.0314,
        .C_prop = 1.0,
        .k_motor = 20.0,
        .k_Tp = 0.0,
        .k_Omega = 0.0,
    };

    return FixedWing(m, Vector4d(Jx, Jy, Jz, Jxz), p);
}

FixedWing generate_sport_drone() {
    // a ~1.4 m span aerobatic plane with a conventional tail, light and stiff in roll
    const double m = 1.8;
    const double Jx = 0.06, Jy = 0.09, Jz = 0.14, Jxz = 0.002;

    const FixedWingParams p{
        .S = 0.35,
        .b = 1.4,
        .c = 0.25,
        .e = 0.85,

        .CL0 = 0.15,
        .CL_alpha = 4.5,
        .CL_q = 5.0,
        .CL_de = 0.4,
        .CD0 = 0.03,
        .CD_alpha = 0.2,
        .CD_q = 0.0,
        .CD_de = 0.02,
        .CD_p = 0.03,
        .Cm0 = 0.0,
        .Cm_alpha = -0.6,
        .Cm_q = -8.0,
        .Cm_de = -0.8,

        .M = 50.0,
        .alpha0 = 0.4712,

        .CY0 = 0.0,
        .CY_beta = -0.5,
        .CY_p = 0.0,
        .CY_r = 0.0,
        .CY_da = 0.0,
        .CY_dr = 0.15,
        .Cl0 = 0.0,
        .Cl_beta = -0.08,
        .Cl_p = -0.5,
        .Cl_r = 0.08,
        .Cl_da = 0.35,
        .Cl_dr = 0.01,
        .Cn0 = 0.0,
        .Cn_beta = 0.12,
        .Cn_p = -0.02,
        .Cn_r = -0.2,
        .Cn_da = 0.01,
        .Cn_dr = -0.08,

        .S_prop = 0.05,
        .C_prop = 1.0,
        .k_motor = 35.0,
        .k_Tp = 0.0,
        .k_Omega = 0.0,
    };

    return FixedWing(m, Vector4d(Jx, Jy, Jz, Jxz), p);
}

Spacecraft generate_cubesat() {
    Spacecraft s;
    SpacecraftParams &p = s.params;

    // empty box with dimensions
    // 0.10 x 0.226 x 0.366 m,
    // m/12 (b^2 + c^2) on the diagonal,
    // plus a bit of inertia on the diagonals (just to see something not symmetric)
    p.m_empty = 11.0;
    p.J_empty << 0.170, -0.001, 0.002, //
        -0.001, 0.130, -0.0005,        //
        0.002, -0.0005, 0.056;

    // tank is on the z axis below the empty cm, so the cm will stay on the thrust line as it gets
    // emptied
    p.r_tank = Vector3d(0.0, 0.0, -0.08);
    p.R_tank = 0.04;

    // small cubesat wheels (quick search says ~10000 rpm at l_max is reasonable)
    p.wheels.A = Matrix3d::Identity();
    p.wheels.J_w = 2.0e-5;
    p.wheels.tau_max = 2.0e-3;
    p.wheels.h_max = 0.02;

    // one main thruster on the -z face pushing along +z
    p.thrusters.push_back(
        {.pos = Vector3d(0.0, 0.0, -0.183), .dir = Vector3d::UnitZ(), .T_max = 1.0, .Isp = 220.0});
    s.controls.throttle.assign(p.thrusters.size(), 0.0);

    const auto [r, v] = circular_orbit(500.0e3, 51.6 * std::numbers::pi / 180.0);
    s.state.pos() = r;
    s.state.vel() = v;
    s.state.mass() = p.m_empty + 1.0; // ssuming a full tank

    return s;
}

namespace {

// where an actuator says the surface is (ideal means no actuator, so the command goes through)
double actual(const Actuator &a, double cmd) { return a.ideal() ? cmd : a.pos; }

} // namespace

void Rocket::step_actuators(double dt) {
    TvcInputs in;
    in.delta_cmd << commands.gimbal_p, commands.gimbal_y;
    in.load = tvc_load;
    in.thrust = tvc_thrust;
    const Eigen::Vector2d rate_before = tvc_state.segment<2>(2);
    tvc_state = tvc_step(tvc, tvc_state, in, dt);
    tvc_accel = dt > 0.0 ? Eigen::Vector2d((tvc_state.segment<2>(2) - rate_before) / dt)
                         : Eigen::Vector2d::Zero();

    for (size_t k = 0; k < fins.size() && k < commands.fins.size(); k++) {
        fins[k].step(commands.grid_fins_deployed ? commands.fins[k] : 0.0, dt);
    }
}

RocketControls Rocket::controls() const {
    RocketControls u = commands;
    // where the TVC actually has the engine, and how it's moving (for the swinging engine's
    // dynamics in the plant)
    u.gimbal_p = tvc_state(0);
    u.gimbal_y = tvc_state(1);
    u.gimbal_p_dot = tvc_state(2);
    u.gimbal_y_dot = tvc_state(3);
    u.gimbal_p_ddot = tvc_accel.x();
    u.gimbal_y_ddot = tvc_accel.y();

    for (size_t k = 0; k < fins.size() && k < u.fins.size(); k++) {
        u.fins[k] = actual(fins[k], commands.grid_fins_deployed ? commands.fins[k] : 0.0);
    }
    return u;
}

Rocket generate_example_rocket() {
    constexpr double PI = std::numbers::pi;
    constexpr double DEG_ = PI / 180.0;

    Rocket r;
    RocketParams &p = r.params;

    // --- geometry ---
    // ICLR's Ember, as built in OpenRocket (data/openrocket/ember.ork): 7.5 m long, 0.4 m diameter,
    // 1.5 m tangent ogive nose on a 6 m tube with a flat base. the body origin is the dry cm, which
    // OpenRocket puts 4.974 m back from the tip, so x = 4.974 - X for any OpenRocket station X
    const double X_cg = 4.974;
    const double x_tip = X_cg, x_tail = X_cg - 7.5, L_nose = 1.5;
    p.d = 0.4;
    p.A_base = p.A_ref();

    // --- mass ---
    // the dry rocket from OpenRocket's component analysis (shells plus the eight internal masses),
    // with the inertias worked out from those components. TODO: placeholder masses until ICLR's
    p.m_dry = 161.5;
    p.x_dry = 0.0;
    p.Ix_dry = 3.07;
    p.It_dry = 572.0;

    // 60 kg of propellant, lox and ethanol at a mixture ratio of 1.3, in the two tanks of the
    // OpenRocket model (lox above, ethanol below, each column sitting on its tank's aft dome)
    p.propellant.push_back({.kind = PropellantStore::Kind::Liquid,
                            .share = 1.3 / 2.3,
                            .x = X_cg - 4.9,
                            .R = 0.18,
                            .rho = 1141.0,
                            .L = 0.0});
    p.propellant.push_back({.kind = PropellantStore::Kind::Liquid,
                            .share = 1.0 / 2.3,
                            .x = X_cg - 6.0,
                            .R = 0.18,
                            .rho = 789.0,
                            .L = 0.0});

    // --- engine ---
    // sized like the fake OpenRocket motor (data/openrocket/ember_fake.eng): 5 kN at sea level,
    // c_eff 2162 m/s, so with a 0.1 m exit that's 5.8 kN vacuum at 256 s. gimballed 0.5 m above the
    // base with the exit flush with it. TODO: Ember's real engine
    RocketEngine &e = p.engine;
    e.kind = RocketEngine::Kind::Throttleable;
    e.A_e = PI * 0.05 * 0.05;
    e.T_vac_max = 5800.0;
    e.Isp_vac = 256.0;
    e.throttle_min = 0.4;
    e.x_G = x_tail + 0.5;
    e.L_n = 0.5;
    p.A_motor = e.A_e;

    e.m_n = 6.0;
    e.l_n = 0.2;
    e.I_na = 0.5 * e.m_n * 0.06 * 0.06;
    e.I_nt = e.m_n * 0.5 * 0.5 / 12.0;

    // --- aero ---
    // the nose is the only barrowman component (CN_alpha = 2, cp 0.696 m back from the tip, both
    // straight from OpenRocket, and constant with mach), the tube adds nothing
    p.components.push_back({.CN_alpha = 2.0, .x_cp = x_tip - 0.696});

    // body outline for the crossflow integral, the nose follows the tangent ogive
    const int n_slices = 75;
    const double dx = (x_tip - x_tail) / n_slices;
    const double R = 0.5 * p.d;
    const double rho_ogive = (R * R + L_nose * L_nose) / (2.0 * R);
    for (int i = 0; i < n_slices; i++) {
        const double x = x_tail + (i + 0.5) * dx;
        const double s = x_tip - x; // distance back from the tip
        const double r_x =
            s >= L_nose ? R
                        : std::sqrt(rho_ogive * rho_ogive - (L_nose - s) * (L_nose - s)) + R -
                              rho_ogive;
        p.slices.push_back({.x = x, .d = 2.0 * r_x, .dx = dx});
    }

    // zero alpha drag against mach, from OpenRocket's component analysis sweep
    // (data/openrocket/ember_cd_*.csv), base drag is left to the model's own formula
    const std::vector<double> mach = {
        0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225,
        0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425,
        0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625,
        0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825,
        0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1, 1.025,
        1.05, 1.075, 1.1, 1.125, 1.15, 1.175, 1.2, 1.225,
        1.25, 1.275, 1.3, 1.325, 1.35, 1.375, 1.4, 1.425,
        1.45, 1.475, 1.5, 1.525, 1.55, 1.575, 1.6, 1.625,
        1.65, 1.675, 1.7, 1.725, 1.75, 1.775, 1.8, 1.825,
        1.85, 1.875, 1.9, 1.925, 1.95, 1.975, 2};
    p.drag.friction = Table1D(mach, {
        0.217, 0.217, 0.217, 0.217, 0.217, 0.217, 0.217, 0.216,
        0.216, 0.216, 0.215, 0.215, 0.215, 0.214, 0.214, 0.213,
        0.213, 0.213, 0.212, 0.211, 0.211, 0.21, 0.21, 0.209,
        0.208, 0.208, 0.207, 0.206, 0.205, 0.204, 0.203, 0.203,
        0.202, 0.201, 0.2, 0.197, 0.194, 0.192, 0.189, 0.187,
        0.184, 0.181, 0.179, 0.177, 0.176, 0.174, 0.173, 0.171,
        0.17, 0.168, 0.167, 0.165, 0.164, 0.162, 0.161, 0.159,
        0.158, 0.156, 0.155, 0.153, 0.152, 0.15, 0.149, 0.147,
        0.146, 0.144, 0.143, 0.142, 0.14, 0.139, 0.137, 0.136,
        0.135, 0.133, 0.132, 0.13, 0.129, 0.128, 0.126});
    p.drag.pressure = Table1D(mach, {
        0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
        0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
        0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
        0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
        0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.036,
        0.058, 0.071, 0.075, 0.072, 0.063, 0.052, 0.038, 0.024,
        0.013, 0.004, 0.0007921, 0.0007574, 0.0007261, 0.0006978, 0.0006721, 0.0006489,
        0.0006276, 0.0006078, 0.0005895, 0.0005726, 0.0005568, 0.000542, 0.000528, 0.000515,
        0.0005027, 0.0004911, 0.0004801, 0.0004697, 0.0004597, 0.0004503, 0.0004413, 0.0004327,
        0.0004245, 0.0004166, 0.0004091, 0.0004018, 0.0003949, 0.0003882, 0.0003818});

    // four grid fins 3.3 m back from the tip (a bit above the middle), in a + (one along each body
    // axis), stowed until the descent (commands.grid_fins_deployed)
    p.fins.n = 4;
    p.fins.x = X_cg - 3.3;
    p.fins.r = R + 0.1;
    p.fins.phi0 = 0.0;
    p.fins.S = 0.2 * 0.2;
    p.fins.CN_alpha = 2.0; // TODO: placeholder
    p.fins.CD0 = 1.2;      // TODO: placeholder

    // --- actuators ---
    // the hydraulic TVC swinging the engine (placeholder geometry and hydraulics from the Ember
    // TVC actuator notes, the nozzle is the engine above)
    r.tvc = tvc_params_for(e);
    r.tvc_state = tvc_rest();

    r.fins.resize(static_cast<size_t>(p.fins.n));
    for (Actuator &f : r.fins) {
        f.kind = Actuator::Kind::SecondOrder;
        f.wn = 30.0;
        f.zeta = 0.7;
        f.rate_max = 150.0 * DEG_;
        f.pos_max = 20.0 * DEG_;
    }

    r.commands.fins.assign(r.fins.size(), 0.0);

    // --- state ---
    // on the pad, nose up (a pitch of +90 deg takes body x to -z, i.e. up), full
    r.state.pos().setZero();
    r.state.vel().setZero();
    r.state.quat() << std::cos(PI / 4.0), 0.0, std::sin(PI / 4.0), 0.0;
    r.state.rates().setZero();
    r.state.mass() = p.m_dry + 60.0;

    return r;
}

} // namespace BSIM
