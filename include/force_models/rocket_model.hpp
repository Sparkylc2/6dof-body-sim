#pragma once
#include "body_state.hpp"
#include "force_models.hpp"
#include "tables.hpp"
#include <Eigen/Dense>
#include <limits>
#include <numbers>
#include <vector>

// --- rocket force model ---
// an axisymmetric rocket with one engine on a two axis gimbal, grid fins
// the aero follows the OpenRocket technical documentation (v13.05) and its source code
// + the variable mass terms left out (jet damping, J dot, coriolis from exhaust) are included
//
// frames (same as the rest of BSIM):
//      inertial is NED, body x points out of the nose, y right, z towards the belly
//      the body origin is a fixed station (the dry cm), and every station x is measured FORWARD
//      along body x from it (openrocket measures X aft from the nose tip, so x = x_tip - X)
//
// the nozzle is treated as massless here, so the thrust just points wherever the gimbal is
namespace BSIM {
using Eigen::Matrix3d;
using Eigen::Vector3d;

constexpr double EARTH_RADIUS_MEAN = 6.371e6; // [m], for gravity falling off with altitude
constexpr double NEVER = std::numeric_limits<double>::infinity(); // an ignition time not reached

// --- propellant ---
// where the propellant sits, each store holds a fixed share of whatever propellant is left
//      liquid settles aft dome (thrust pushes it back) and fills forward as a column
//      solid is a case bonded grain burning outwards from the bore
struct PropellantStore {
    enum class Kind { Liquid, Solid };
    Kind kind = Kind::Liquid;

    double share = 1.0; // fraction of the propellant
    double x = 0.0;   // liquid: station of the aft dome (bottom of the column), solid: grain centre
    double R = 0.0;   // liquid: tank radius, solid: outer radius [m]
    double rho = 0.0; // propellant density [kg/m^3]
    double L = 0.0;   // solid only: grain length [m]
};

// --- engine ---
//  ThrustCurve:  T(t since ignition) from a static fire (sea level), mass flow from c_eff
//  Throttleable: T = throttle T_vac_max - p_a A_e, mass flow from the vacuum Isp
struct RocketEngine {
    enum class Kind { ThrustCurve, Throttleable };
    Kind kind = Kind::Throttleable;

    // thrust curve engines
    Table1D thrust_curve;     // T [N] vs time since ignition [s] (at p_test)
    double c_eff = 0.0;       // effective exhaust velocity, total impulse / propellant mass [m/s]
    double p_test = 101325.0; // ambient pressure from measurement time [Pa]

    // throttleable engines
    double T_vac_max = 0.0;    // full throttle vacuum thrust [N]
    double Isp_vac = 0.0;      // vacuum specific impulse [s]
    double throttle_min = 0.0; // lowest throttle it can run at (below it, it's off)

    // both
    double A_e = 0.0; // nozzle exit area [m^2]
    double x_G = 0.0; // gimbal pivot station [m]
    double L_n = 0.0; // pivot to nozzle exit plane, along the nozzle [m]

    // the part that swings on the gimbal (N in the TVC notes, section 09), see tvc_model.hpp
    // it's already counted in the dry mass at delta = 0, these are only for what changes when it
    // swings. m_n = 0 is the massless nozzle
    double m_n = 0.0;  // its mass [kg]
    double l_n = 0.0;  // its cm behind the pivot, along the nozzle [m]
    double I_na = 0.0; // its inertia about its own axis, about its cm [kg m^2]
    double I_nt = 0.0; // its transverse inertia about its cm [kg m^2]
};

// --- aero ---
// a barrowman body component (nose, shoulder, boattail, or a fin set) at its own station
//      CN_alpha = 2 (A(l) - A(0)) / A_ref for a body component (OpenRocket eq 3.19), per rad
struct BodyComponent {
    double CN_alpha = 0.0; // normal force slope on A_ref [1/rad]
    double x_cp = 0.0;     // station it acts at [m]
};

// one slice of the body for the crossflow strip integral
struct BodySlice {
    double x = 0.0;  // station of the slice middle [m]
    double d = 0.0;  // diameter there [m]
    double dx = 0.0; // slice length [m]
};

// zero angle of attack drag, each on A_ref, against mach
// (an empty table counts as zero, the base drag is worked out from the base area instead)
// TODO: these get filled from proper curves later
struct DragCurves {
    Table1D friction;  // skin friction (OpenRocket 3.4.2)
    Table1D pressure;  // nose / fin / body pressure drag (3.4.3, 3.4.4)
    Table1D parasitic; // launch lugs, rail buttons etc (3.4.6)
};

// grid fins, n of them evenly round the body at one station
// fin k sits at roll angle phi_k = phi0 + 2 pi k / n (phi = 0 is along +y, pi/2 along +z)
// a positive deflection of fin k pushes it along -n_k, where n_k = x x e_r is the fin's
// tangential normal
//
// force flips sign when flow comes from behind
// only act once RocketControls::grid_fins_deployed is set, before that they're stowed flat
// against the body and add nothing
struct GridFins {
    int n = 0;             // number of fins (0 = none)
    double x = 0.0;        // station of the fin centres [m]
    double r = 0.0;        // radius of the fin centre from the body axis [m]
    double phi0 = 0.0;     // roll angle of fin 0 [rad]
    double S = 0.0;        // reference area of one fin [m^2]
    double CN_alpha = 0.0; // normal force slope on S [1/rad]
    double CD0 = 0.0;      // drag on S
};

struct RocketParams {
    // --- mass ---
    double m_dry = 0.0;  // everything but the propellant [kg]
    double x_dry = 0.0;  // its cm (0 if the body origin is the dry cm) [m]
    double Ix_dry = 0.0; // roll inertia about body axis [kg m^2]
    double It_dry = 0.0; // transverse inertia about dry cm [kg m^2]
    std::vector<PropellantStore> propellant;

    RocketEngine engine;

    // --- geometry ---
    double d = 0.0;       // reference (body) diameter [m]
    double A_base = 0.0;  // area of the flat base [m^2]
    double A_motor = 0.0; // base area of the plume while burning (gets taken out of base drag)

    // --- aero ---
    std::vector<BodyComponent> components; // barrowman normal force pieces
    std::vector<BodySlice> slices;         // body outline for the crossflow strip integral
    double C_dc = 1.1;                     // crossflow drag of a cylinder (OpenRocket's K)
    DragCurves drag;
    double CD_flat = 1.17; // flat plate drag, for the base when flying tail first (TODO)
    GridFins fins;

    double A_ref() const { return 0.25 * std::numbers::pi * d * d; }
};

// what the actuators are actually doing (the physics sees these, not the commands)
struct RocketControls {
    double gimbal_p = 0.0;           // pitch gimbal [rad], + gives a nose up (+q) moment
    double gimbal_y = 0.0;           // yaw gimbal [rad], + gives a nose right (+r) moment
    double gimbal_p_dot = 0.0;       // how fast the gimbal is moving [rad/s] (from the servos,
    double gimbal_y_dot = 0.0;       //  only the gimbal dynamics use these and the accelerations)
    double gimbal_p_ddot = 0.0;      // [rad/s^2]
    double gimbal_y_ddot = 0.0;
    double throttle = 0.0;           // [0, 1] (throttleable engines)
    double t_ignition = NEVER;       // sim time the motor was lit (thrust curve engines)
    std::vector<double> fins;        // grid fin deflections [rad], one per fin
    bool grid_fins_deployed = false; // stowed (no aero) until this is set
};

// which pieces are switched on
struct RocketEnvironment {
    bool gravity = true;
    bool aero = true;
    bool jet_damping = true; // -m_dot_e rho_e x (w x rho_e)
    bool J_dot = true;       // -J_dot w
    bool coriolis = true;    // -2 m_dot_e w x rho_e
};

// the air the rocket is in
struct RocketAir {
    double rho = 1.225;                 // density [kg/m^3]
    double a = 340.29399;               // speed of sound [m/s]
    double p = 101325.0;                // ambient pressure [Pa]
    Vector3d wind_i = Vector3d::Zero(); // steady wind (NED)
    Vector3d gust_b = Vector3d::Zero(); // turbulence (body)
};

// --- mass properties ---
struct RocketMassProperties {
    double m = 0.0;
    double x_cm = 0.0; // station of the current cm
    double Ix = 0.0;   // about the body axis
    double It = 0.0;   // transverse, about the current cm

    Vector3d r_cm() const { return Vector3d(x_cm, 0.0, 0.0); }
    Matrix3d J() const { return Vector3d(Ix, It, It).asDiagonal(); }
};

// mass, cm and inertia for a total mass m (propellant = m - m_dry, held at >= 0)
RocketMassProperties rocket_mass_properties(const RocketParams &p, double m);

// J dot from a central difference of the mass properties in m, times m_dot
Matrix3d rocket_J_dot(const RocketParams &p, double m, double m_dot);

// --- engine ---
struct EngineOutput {
    double T = 0.0;                    // thrust [N]
    double m_dot_e = 0.0;              // exhaust mass flow [kg/s] (m_dot = -m_dot_e)
    Vector3d n_T = Vector3d::UnitX();  // direction it pushes the rocket (body)
    Vector3d rho_e = Vector3d::Zero(); // nozzle exit from the cm (body) [m]
    Wrench wrench;                     // thrust force and moment about the cm (body)
};

// unit thrust direction for the gimbal angles, R_y(-dp) R_z(-dy) x
Vector3d thrust_direction(double gimbal_p, double gimbal_y);

// thrust, mass flow and the thrust wrench about the cm at x_cm
// (off once the propellant is gone, or before ignition / below minimum throttle)
EngineOutput rocket_engine(const RocketParams &p, const RocketControls &u, double t, double m,
                           double x_cm, double p_ambient);

// --- air data ---
struct RocketAirData {
    Vector3d v_r = Vector3d::Zero(); // air relative velocity, body [u_r v_r w_r]
    double V = 0.0;                  // airspeed
    double v_perp = 0.0;             // crossflow speed sqrt(v_r^2 + w_r^2)
    double alpha_t = 0.0;            // total angle of attack, in [0, pi]
    double mach = 0.0;
    double qbar = 0.0;
    double rho = 0.0;
};

RocketAirData rocket_air_data(const RocketState &x, const RocketAir &air);

// --- aero ---
// OpenRocket's axial drag multiplier on [0, pi/2] (BarrowmanCalculator.calculateAxialCD)
//      1 + 0.3 (3s^2 - 2s^3) up to 17 deg, then 1.3 (1-u)^3 (1+3u) down to 0 at 90 deg
double axial_drag_multiplier(double alpha);

// OpenRocket's base drag coefficient (eq 3.94), on the base area
double base_drag_coefficient(double mach);

// zero alpha drag on A_ref, the curves plus base drag on the uncovered base
double rocket_CD0(const RocketParams &p, double mach, bool burning);

// axial force coefficient on A_ref (positive pushes the rocket backwards, along -x)
//      nose first: CD0 * multiplier (openrocket)
//      tail first (alpha_t > 90 deg): flat plate drag on the base, pushing along +x
//      TODO: proper tail first data
double rocket_CA(const RocketParams &p, double mach, bool burning, double alpha_t);

// the pieces, each a wrench about the cm at x_cm (body frame)
Wrench rocket_axial_force(const RocketParams &p, const RocketAirData &ad, bool burning);
Wrench rocket_barrowman_normal(const RocketParams &p, const RocketAirData &ad, double x_cm);
// crossflow drag of every slice with its own local crossflow (v_r, w_r) + w x xi, which gives
// the body lift and the (quadratic) pitch/yaw damping together
Wrench rocket_crossflow(const RocketParams &p, const RocketAirData &ad, const Vector3d &w,
                        double x_cm);
Wrench rocket_grid_fins(const RocketParams &p, const RocketControls &u, const RocketAirData &ad,
                        const Vector3d &w, double x_cm);

// everything aerodynamic
Wrench rocket_aero(const RocketParams &p, const RocketControls &u, const RocketAirData &ad,
                   const Vector3d &w, double x_cm, bool burning);

// --- gravity ---
// falls off w alt h = -pd
Vector3d rocket_gravity(double h);

// --- gimbal dynamics ---
// every moment the swinging gimbal mass adds, from rocket_gimbal_dynamics (tvc_model.hpp)
using Matrix32d = Eigen::Matrix<double, 3, 2>;
struct GimbalDynamics {
    Matrix3d dJ = Matrix3d::Zero();      // J(m, delta) - J(m, 0)
    Matrix32d B = Matrix32d::Zero();
    Vector3d H_rel = Vector3d::Zero();   // B delta_dot
    Vector3d M_twd = Vector3d::Zero();   // -B delta_ddot
    Vector3d M_cor = Vector3d::Zero();   // -w x B delta_dot - (dJ/ddelta delta_dot) w
    Vector3d M_geom = Vector3d::Zero();  // -B_dot delta_dot - (dJ/dx_cm x_cm_dot) w
    Vector3d M_exit = Vector3d::Zero();  // -m_dot_e rho_e x rho_e_dot

    Vector3d M() const { return M_twd + M_cor + M_geom + M_exit; }
};

// --- everything ---
// every term in the equations of motion (body frame unless it says otherwise), for
// logging/accelerometer model
struct RocketTerms {
    RocketMassProperties mp;
    Matrix3d J_dot = Matrix3d::Zero();
    EngineOutput engine;
    RocketAirData air;
    Wrench aero;
    Vector3d F_cor = Vector3d::Zero();   // -2 m_dot_e w x rho_e
    Vector3d M_jet = Vector3d::Zero();   // -m_dot_e rho_e x (w x rho_e)
    Vector3d M_J_dot = Vector3d::Zero(); // -J_dot w
    Vector3d g_i = Vector3d::Zero();     // gravity (NED)
    GimbalDynamics gimbal;               // the swinging gimbal mass (zero with a massless nozzle)
    Eigen::Vector2d hinge = Eigen::Vector2d::Zero(); // servo hinge torques [pitch, yaw] [N m]
    Eigen::Vector2d nozzle_load = Eigen::Vector2d::Zero(); // the vehicle's load on the nozzle,
                                                           // the TVC's TvcInputs::load [N m]
                                                     // (filled in by rocket_derivatives)

    // total non gravitational force, what an accelerometer at the cm reads times m
    Vector3d F_body() const { return engine.wrench.F + aero.F + F_cor; }
    // total moment about the cm
    Vector3d M_body() const { return engine.wrench.M + aero.M + M_jet + M_J_dot + gimbal.M(); }
    // inertia about the cm with the gimbal where it is
    Matrix3d J() const { return mp.J() + gimbal.dJ; }
};

RocketTerms rocket_terms(const RocketParams &p, const RocketEnvironment &env,
                         const RocketControls &u, const RocketAir &air, const RocketState &x,
                         double t);

} // namespace BSIM
