#include "body.hpp"
#include <algorithm>
#include <raylib.h>

namespace uav {

// --- drone ---
Drone::Drone() {
    m = 1.56;

    const double Ix = 0.1147, Iy = 0.0576, Iz = 0.1712, Ixz = 0.0015;

    // following the symmetry assumption that airplanes use
    I << Ix, Iy, Iz, Ixz;
}

void Drone::read_keyboard() {
    input.F = keyboard_force();
    input.M = keyboard_moment();
}

Vector3d Drone::keyboard_force() const {
    Vector3d f = Vector3d::Zero();

    double mult = 2.0 * m * 9.81;
    if (IsKeyDown(KEY_RIGHT_SHIFT)) {
        mult *= 1.2;
    }

    if (IsKeyDown(KEY_SPACE)) {
        f += Vector3d(0.0, 0.0, -mult);
    }
    if (IsKeyDown(KEY_LEFT_SHIFT)) {
        f += Vector3d(0.0, 0.0, mult);
    }

    return f;
}

Vector3d Drone::keyboard_moment() const {
    Vector3d T = Vector3d::Zero();

    double mult = 0.5;

    if (IsKeyDown(KEY_RIGHT_SHIFT)) {
        mult *= 1.2;
    }
    if (IsKeyDown(KEY_W)) {
        T.y() -= mult;
    }
    if (IsKeyDown(KEY_S)) {
        T.y() += mult;
    }
    if (IsKeyDown(KEY_A)) {
        T.x() -= mult;
    }
    if (IsKeyDown(KEY_D)) {
        T.x() += mult;
    }
    if (IsKeyDown(KEY_Q)) {
        T.z() -= mult;
    }
    if (IsKeyDown(KEY_E)) {
        T.z() += mult;
    }

    return T;
}

Wrench Drone::wrench(const AircraftState &x, const AirState &) const {
    Wrench w = gravity(m, x.quat());
    w += input;
    return w;
}

// --- fixed wing ---
FixedWing::FixedWing(double m, const Vector4d &I, const FixedWingParams &params) : params(params) {
    this->m = m;
    this->I = I;
}

void FixedWing::read_keyboard(double dt) {
    const double elevator = 0.1;      // rad on top of trim
    const double aileron = 0.1;       // rad on top of trim
    const double rudder = 0.15;       // rad on top of trim
    const double throttle_rate = 0.5; // per second

    // positive elevator pitches the nose down (Cm_de < 0)
    controls.de = trim.de;
    if (IsKeyDown(KEY_W)) {
        controls.de += elevator;
    }
    if (IsKeyDown(KEY_S)) {
        controls.de -= elevator;
    }

    // positive aileron rolls right (Cl_da > 0)
    controls.da = trim.da;
    if (IsKeyDown(KEY_D)) {
        controls.da += aileron;
    }
    if (IsKeyDown(KEY_A)) {
        controls.da -= aileron;
    }

    // positive rudder yaws the nose left (Cn_dr < 0)
    controls.dr = trim.dr;
    if (IsKeyDown(KEY_Q)) {
        controls.dr += rudder;
    }
    if (IsKeyDown(KEY_E)) {
        controls.dr -= rudder;
    }

    // throttle stays where it's left
    if (IsKeyDown(KEY_UP)) {
        controls.dt += throttle_rate * dt;
    }
    if (IsKeyDown(KEY_DOWN)) {
        controls.dt -= throttle_rate * dt;
    }
    controls.dt = std::clamp(controls.dt, 0.0, 1.0);
}

Wrench FixedWing::wrench(const AircraftState &x, const AirState &air) const {

    const AirData air_info = air_data(x.vel(), x.quat(), air);

    Wrench w = gravity(m, x.quat());
    w += longitudinal_aero(params, air_info, x.rates().y(), controls.de);
    w += lateral_aero(params, air_info, x.rates(), controls.da, controls.dr);
    w += propulsion(params, air_info, controls.dt);
    return w;
}

// factories for the fixed wing
FixedWing generate_zagi_flying_wing() {
    const double m = 1.56;
    const double Ix = 0.1147, Iy = 0.0576, Iz = 0.1712, Ixz = 0.0015;

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

    return FixedWing(m, Vector4d(Ix, Iy, Iz, Ixz), p);
}

FixedWing generate_sport_drone() {
    // a ~1.4 m span aerobatic plane with a conventional tail, light and stiff in roll
    const double m = 1.8;
    const double Ix = 0.06, Iy = 0.09, Iz = 0.14, Ixz = 0.002;

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

    return FixedWing(m, Vector4d(Ix, Iy, Iz, Ixz), p);
}

} // namespace uav
