#include "force_models/fixed_wing_model.hpp"
#include <numbers>

namespace BSIM {

double stall_blend(const FixedWingParams &p, double alpha) {
    // alpha is in [-pi, pi] (from atan2) so with M around 50 these stay well inside a double
    const double a = std::exp(-p.M * (alpha - p.alpha0));
    const double b = std::exp(p.M * (alpha + p.alpha0));
    return (1.0 + a + b) / ((1.0 + a) * (1.0 + b));
}

double lift_coefficient(const FixedWingParams &p, double alpha) {
    const double sigma = stall_blend(p, alpha);
    const double linear = p.CL0 + p.CL_alpha * alpha;
    const double sign = alpha >= 0.0 ? 1.0 : -1.0;
    const double s = std::sin(alpha), c = std::cos(alpha);

    // linear lift before stall, flat plate lift after
    return (1.0 - sigma) * linear + sigma * (2.0 * sign * s * s * c);
}

double drag_coefficient(const FixedWingParams &p, double alpha) {
    const double AR = p.b * p.b / p.S; // aspect ratio
    const double sigma = stall_blend(p, alpha);
    const double CL_linear = p.CL0 + p.CL_alpha * alpha;
    const double s = std::sin(alpha);

    // parasitic + induced
    const double polar = p.CD_p + CL_linear * CL_linear / (std::numbers::pi * p.e * AR);
    // the flat plate normal force is 2sin^2(alpha) (same plate lift is from),
    // so the component along the freestream is the drag
    const double flat_plate = 2.0 * s * s * std::abs(s);

    // polar before stall, then flat plate drag after
    return (1.0 - sigma) * polar + sigma * flat_plate;
}

Wrench longitudinal_aero(const FixedWingParams &p, const AirData &air, double q, double de) {
    const double ca = std::cos(air.alpha), sa = std::sin(air.alpha);

    const double CL = lift_coefficient(p, air.alpha);
    const double CD = drag_coefficient(p, air.alpha);

    // lift and drag rotated from the stability frame into body x and z
    const double CX = -CD * ca + CL * sa;
    const double CX_q = -p.CD_q * ca + p.CL_q * sa;
    const double CX_de = -p.CD_de * ca + p.CL_de * sa;
    const double CZ = -CD * sa - CL * ca;
    const double CZ_q = -p.CD_q * sa - p.CL_q * ca;
    const double CZ_de = -p.CD_de * sa - p.CL_de * ca;

    // dynamic pressure times wing area
    const double qbar_S = 0.5 * air.rho * air.Va * air.Va * p.S;
    // qbar_S * c / (2 Va) multiplied out so the Va cancels and it's fine at rest
    const double qbar_S_c_2Va = 0.25 * air.rho * air.Va * p.S * p.c;

    Wrench w;
    w.F.x() = qbar_S * (CX + CX_de * de) + qbar_S_c_2Va * CX_q * q;
    w.F.z() = qbar_S * (CZ + CZ_de * de) + qbar_S_c_2Va * CZ_q * q;
    w.M.y() = p.c * (qbar_S * (p.Cm0 + p.Cm_alpha * air.alpha + p.Cm_de * de) +
                     qbar_S_c_2Va * p.Cm_q * q);
    return w;
}

Wrench lateral_aero(const FixedWingParams &p, const AirData &air, const Vector3d &rates, double da,
                    double dr) {
    const double roll_rate = rates.x(), yaw_rate = rates.z();
    const double beta = air.beta;

    const double qbar_S = 0.5 * air.rho * air.Va * air.Va * p.S;
    const double qbar_S_b_2Va = 0.25 * air.rho * air.Va * p.S * p.b;

    Wrench w;
    // these all just follow from a taylor series expansion
    w.F.y() = qbar_S * (p.CY0 + p.CY_beta * beta + p.CY_da * da + p.CY_dr * dr) +
              qbar_S_b_2Va * (p.CY_p * roll_rate + p.CY_r * yaw_rate);
    w.M.x() = p.b * (qbar_S * (p.Cl0 + p.Cl_beta * beta + p.Cl_da * da + p.Cl_dr * dr) +
                     qbar_S_b_2Va * (p.Cl_p * roll_rate + p.Cl_r * yaw_rate));
    w.M.z() = p.b * (qbar_S * (p.Cn0 + p.Cn_beta * beta + p.Cn_da * da + p.Cn_dr * dr) +
                     qbar_S_b_2Va * (p.Cn_p * roll_rate + p.Cn_r * yaw_rate));
    return w;
}

Wrench propulsion(const FixedWingParams &p, const AirData &air, double dt) {
    const double prop_speed = p.k_motor * dt;
    const double prop_rate = p.k_Omega * dt;

    Wrench w;
    // goes negative (windmilling drag) once Va is faster than the air leaving the prop
    w.F.x() = 0.5 * air.rho * p.S_prop * p.C_prop * (prop_speed * prop_speed - air.Va * air.Va);
    w.M.x() = -p.k_Tp * prop_rate * prop_rate;
    return w;
}
} // namespace BSIM
