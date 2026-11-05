#pragma once
#include "force_models.hpp"
#include <Eigen/Dense>

namespace BSIM {
// gravity on a body of mass m with attitude e (quaternion)
// acts at the center of mass so there's no moment
Wrench gravity(double m, const Vector4d &e);

// --- fixed wing (chapter 4) ---
// aero/propulsive coefficients for a fixed wing
// everything is SI and per radian
struct FixedWingParams {
    // geometry (air density comes from the AirData now)
    double S = 0.0; // wing area
    double b = 0.0; // wing span
    double c = 0.0; // mean chord
    double e = 0.0; // oswald efficiency (for the nonlinear drag)

    // longitudinal (lift, drag, pitch moment)
    double CL0 = 0.0, CL_alpha = 0.0, CL_q = 0.0, CL_de = 0.0;
    double CD0 = 0.0, CD_alpha = 0.0, CD_q = 0.0, CD_de = 0.0;
    double CD_p = 0.0; // parasitic drag (for the nonlinear drag)
    double Cm0 = 0.0, Cm_alpha = 0.0, Cm_q = 0.0, Cm_de = 0.0;

    // stall blending between linear and flat plate lift
    double M = 0.0;      // how sharp the transition is
    double alpha0 = 0.0; // stall angle

    // lateral (side force, roll moment, yaw moment)
    double CY0 = 0.0, CY_beta = 0.0, CY_p = 0.0, CY_r = 0.0, CY_da = 0.0, CY_dr = 0.0;
    double Cl0 = 0.0, Cl_beta = 0.0, Cl_p = 0.0, Cl_r = 0.0, Cl_da = 0.0, Cl_dr = 0.0;
    double Cn0 = 0.0, Cn_beta = 0.0, Cn_p = 0.0, Cn_r = 0.0, Cn_da = 0.0, Cn_dr = 0.0;

    // propulsion
    double S_prop = 0.0, C_prop = 0.0, k_motor = 0.0, k_Tp = 0.0, k_Omega = 0.0;
};

// control surface deflections in radians and throttle in [0, 1]
struct FixedWingControls {
    double de = 0.0; // elevator
    double da = 0.0; // aileron
    double dr = 0.0; // rudder
    double dt = 0.0; // throttle
};

// sigmoid blending linear lift into the flat plate lift around the stall angle
// from eqn 4.10
double stall_blend(const FixedWingParams &p, double alpha);

// nonlinear lift and drag coefficients, only the alpha dependent part
// from eqn 4.9
double lift_coefficient(const FixedWingParams &p, double alpha);
// from eqn. 4.11, blended into flat plate drag past stall with the same sigmoid as the lift
double drag_coefficient(const FixedWingParams &p, double alpha);

// lift, drag (rotated into body x and z) and pitch moment, with q is the pitch rate
// from eqn. 4.18 (lift and drag), and eqn. 4.20 (pitching moment)
Wrench longitudinal_aero(const FixedWingParams &p, const AirData &air, double q, double de);

// side force, roll and yaw moment, rates are the body rates [p, q, r]
// from eqn 4.18 (side force), and eqn 4.20 (roll and yaw moments)
Wrench lateral_aero(const FixedWingParams &p, const AirData &air, const Vector3d &rates, double da,
                    double dr);

// propeller thrust along body x and its torque about body x
// this assumes it the thrustline goes through the centroid
// from eqn 4.18 (thrust force), and eqn 4.20 (prop gyro moment)
Wrench propulsion(const FixedWingParams &p, const AirData &air, double dt);
} // namespace BSIM
