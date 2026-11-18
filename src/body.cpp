#include "body.hpp"

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

} // namespace BSIM
