#include "force_models/rocket_model.hpp"
#include "atmosphere.hpp" // ISA_G0
#include "frames.hpp"
#include <algorithm>
#include <cmath>

namespace BSIM {

namespace {

constexpr double PI = std::numbers::pi;
constexpr double ALPHA_BREAKPOINT = 17.0 * PI / 180.0; // openrocket's axial drag breakpoint

// an empty curve is just 0
double curve(const Table1D &t, double x) { return t.x.empty() ? 0.0 : t(x); }

// inertia about its own cm
struct PropMass {
    double m = 0.0, x = 0.0, Ix = 0.0, It = 0.0;
};

PropMass store_mass(const PropellantStore &s, double m_prop) {
    PropMass o;
    o.m = s.share * m_prop;
    o.x = s.x;

    if (o.m <= 0.0) {
        return o;
    }

    if (s.kind == PropellantStore::Kind::Liquid) {
        // column of height h on the aft dome
        const double h = o.m / (s.rho * PI * s.R * s.R);
        o.x = s.x + 0.5 * h;
        o.It = o.m * (3.0 * s.R * s.R + h * h) / 12.0;
    } else {
        // hollow cylinder, the bore grows as it burns
        const double Ri2 = std::max(s.R * s.R - o.m / (s.rho * PI * s.L), 0.0);
        o.Ix = 0.5 * o.m * (s.R * s.R + Ri2);
        o.It = o.m * ((s.R * s.R + Ri2) / 4.0 + s.L * s.L / 12.0);
    }

    return o;
}

} // namespace

// --- mass properties ---
RocketMassProperties rocket_mass_properties(const RocketParams &p, double m) {
    const double m_prop = std::max(m - p.m_dry, 0.0);

    RocketMassProperties mp;
    mp.m = p.m_dry + m_prop;

    // cm first, the transverse inertia needs it
    double mx = p.m_dry * p.x_dry;
    for (const PropellantStore &s : p.propellant) {
        const PropMass sm = store_mass(s, m_prop);
        mx += sm.m * sm.x;
    }

    mp.x_cm = mx / mp.m;

    // parallel axis for everything about the current cm
    mp.Ix = p.Ix_dry;
    mp.It = p.It_dry + p.m_dry * (p.x_dry - mp.x_cm) * (p.x_dry - mp.x_cm);

    for (const PropellantStore &s : p.propellant) {
        const PropMass sm = store_mass(s, m_prop);
        mp.Ix += sm.Ix;
        mp.It += sm.It + sm.m * (sm.x - mp.x_cm) * (sm.x - mp.x_cm);
    }

    return mp;
}

Matrix3d rocket_J_dot(const RocketParams &p, double m, double m_dot) {
    if (m_dot == 0.0) {
        return Matrix3d::Zero();
    }

    // central difference, one sided at empty (no propellant below m_dry to difference into)
    constexpr double EPS = 1e-3; // [kg]
    const double lo = std::max(m - EPS, p.m_dry);
    const double hi = m + EPS;

    const Matrix3d dJ_dm =
        (rocket_mass_properties(p, hi).J() - rocket_mass_properties(p, lo).J()) / (hi - lo);
    return dJ_dm * m_dot;
}

// --- engine ---
Vector3d thrust_direction(double gimbal_p, double gimbal_y) {
    return Vector3d(std::cos(gimbal_p) * std::cos(gimbal_y), -std::sin(gimbal_y),
                    std::sin(gimbal_p) * std::cos(gimbal_y));
}

EngineOutput rocket_engine(const RocketParams &p, const RocketControls &u, double t, double m,
                           double x_cm, double p_ambient) {
    const RocketEngine &e = p.engine;
    const Vector3d r_G(e.x_G, 0.0, 0.0);
    const Vector3d r_cm(x_cm, 0.0, 0.0);

    EngineOutput out;
    out.n_T = thrust_direction(u.gimbal_p, u.gimbal_y);

    out.rho_e = r_G - e.L_n * out.n_T - r_cm;

    // nothing left to burn
    if (m <= p.m_dry) {
        return out;
    }

    switch (e.kind) {
    case RocketEngine::Kind::ThrustCurve: {
        const double t_burn = t - u.t_ignition;
        if (e.thrust_curve.x.empty() || t_burn < 0.0 || t_burn > e.thrust_curve.x.back()) {
            return out;
        }

        // curve is measured at p_test, so it gains (p_test - p_a) A_e going up
        // (mass flow doesn't change with the ambient pressure, only the pressure thrust)
        const double T_curve = std::max(e.thrust_curve(t_burn), 0.0);
        out.T = std::max(T_curve + (e.p_test - p_ambient) * e.A_e, 0.0);
        out.m_dot_e = e.c_eff > 0.0 ? T_curve / e.c_eff : 0.0;
        break;
    }

    case RocketEngine::Kind::Throttleable: {
        const double tau = std::clamp(u.throttle, 0.0, 1.0);
        // below its minimum throttle it can't run, so it's off
        if (tau <= 0.0 || tau < e.throttle_min) {
            return out;
        }

        const double T_vac = tau * e.T_vac_max;
        out.T = std::max(T_vac - p_ambient * e.A_e, 0.0);
        out.m_dot_e = T_vac / (ISA_G0 * e.Isp_vac);
        break;
    }
    }

    out.wrench.F = out.T * out.n_T;
    out.wrench.M = (r_G - r_cm).cross(out.wrench.F);
    return out;
}

// --- air data ---
RocketAirData rocket_air_data(const RocketState &x, const RocketAir &air) {
    RocketAirData ad;
    ad.rho = air.rho;
    ad.v_r = rotation(x.quat()).transpose() * (x.vel() - air.wind_i) - air.gust_b;
    ad.V = ad.v_r.norm();
    ad.v_perp = std::hypot(ad.v_r.y(), ad.v_r.z());
    // atan2 gives 0 at rest, so nothing below needs a special case there
    ad.alpha_t = std::atan2(ad.v_perp, ad.v_r.x());
    ad.mach = air.a > 0.0 ? ad.V / air.a : 0.0;
    ad.qbar = 0.5 * air.rho * ad.V * ad.V;
    return ad;
}

// --- aero ---
double axial_drag_multiplier(double alpha) {
    alpha = std::clamp(alpha, 0.0, 0.5 * PI);

    if (alpha <= ALPHA_BREAKPOINT) {
        const double s = alpha / ALPHA_BREAKPOINT;
        return 1.0 + 0.3 * (3.0 * s * s - 2.0 * s * s * s);
    }

    const double u = (alpha - ALPHA_BREAKPOINT) / (0.5 * PI - ALPHA_BREAKPOINT);
    return 1.3 * (1.0 - u) * (1.0 - u) * (1.0 - u) * (1.0 + 3.0 * u);
}

double base_drag_coefficient(double mach) {
    return mach < 1.0 ? 0.12 + 0.13 * mach * mach : 0.25 / mach;
}

double rocket_CD0(const RocketParams &p, double mach, bool burning) {
    const double A_ref = p.A_ref();
    if (A_ref <= 0.0) {
        return 0.0;
    }

    // plume covers the base so we can get rid of it
    const double A_base_eff = std::max(p.A_base - (burning ? p.A_motor : 0.0), 0.0);

    return curve(p.drag.friction, mach) + curve(p.drag.pressure, mach) +
           curve(p.drag.parasitic, mach) + base_drag_coefficient(mach) * A_base_eff / A_ref;
}

double rocket_CA(const RocketParams &p, double mach, bool burning, double alpha_t) {
    if (alpha_t <= 0.5 * PI) {
        return rocket_CD0(p, mach, burning) * axial_drag_multiplier(alpha_t);
    }

    // TODO: tail first. for now the base is a flat plate facing the flow, its drag projected on
    // the axis (-cos alpha_t goes 0 -> 1 from 90 to 180 deg)
    const double A_ref = p.A_ref();
    if (A_ref <= 0.0) {
        return 0.0;
    }
    return -p.CD_flat * (p.A_base / A_ref) * (-std::cos(alpha_t));
}

Wrench rocket_axial_force(const RocketParams &p, const RocketAirData &ad, bool burning) {
    Wrench w;
    w.F.x() = -ad.qbar * p.A_ref() * rocket_CA(p, ad.mach, burning, ad.alpha_t);
    return w;
}

Wrench rocket_barrowman_normal(const RocketParams &p, const RocketAirData &ad, double x_cm) {
    // qbar A_ref CN_alpha sin(alpha_t) against against da crossflow
    const double k = -0.5 * ad.rho * p.A_ref() * ad.V;

    Wrench w;
    for (const BodyComponent &c : p.components) {
        const double Fy = k * c.CN_alpha * ad.v_r.y();
        const double Fz = k * c.CN_alpha * ad.v_r.z();
        const double Delta = x_cm - c.x_cp; // + if the component is behind cm

        w.F += Vector3d(0.0, Fy, Fz);
        w.M += Vector3d(0.0, Delta * Fz, -Delta * Fy);
    }
    return w;
}

Wrench rocket_crossflow(const RocketParams &p, const RocketAirData &ad, const Vector3d &w,
                        double x_cm) {
    Wrench out;

    for (const BodySlice &s : p.slices) {
        // the curr slice moves at v_r + w x (xi x)
        const double xi = s.x - x_cm;
        const double cy = ad.v_r.y() + w.z() * xi;
        const double cz = ad.v_r.z() - w.y() * xi;
        const double c = std::hypot(cy, cz);

        // cylinder in crossflow, 1/2 rho C_dc d |c| c per unit length against the crossflow
        const double k = -0.5 * ad.rho * p.C_dc * s.d * c * s.dx;
        const double Fy = k * cy;
        const double Fz = k * cz;

        out.F += Vector3d(0.0, Fy, Fz);
        out.M += Vector3d(0.0, -xi * Fz, xi * Fy); // (xi x) x (0, Fy, Fz)
    }

    return out;
}

Wrench rocket_grid_fins(const RocketParams &p, const RocketControls &u, const RocketAirData &ad,
                        const Vector3d &w, double x_cm) {
    const GridFins &g = p.fins;
    Wrench out;
    // stowed fins sit flat against the body, so they do nothing
    if (g.n <= 0 || !u.grid_fins_deployed) {
        return out;
    }

    const Vector3d r_cm(x_cm, 0.0, 0.0);

    for (int k = 0; k < g.n; k++) {
        const double phi = g.phi0 + 2.0 * PI * k / g.n;
        const double c = std::cos(phi), s = std::sin(phi);

        const Vector3d n_k(0.0, -s, c);            // fin normal (tangential), x cross e_r
        const Vector3d r_k(g.x, g.r * c, g.r * s); // fin centre
        const Vector3d v_l = ad.v_r + w.cross(r_k - r_cm);

        const size_t ku = static_cast<size_t>(k);
        const double delta = ku < u.fins.size() ? u.fins[ku] : 0.0;

        // deflecting by delta about the fin's radial axis tilts its normal towards x
        const Vector3d m_k = n_k * std::cos(delta) + Vector3d::UnitX() * std::sin(delta);

        // linear normal force (sin of the local incidence times V^2, with the V's cancelled),
        // the same whichever way the flow comes from, plus drag against the local flow
        const double V_l = v_l.norm();
        const double v_m = v_l.dot(m_k);
        const Vector3d F = -0.5 * ad.rho * g.S * (g.CN_alpha * V_l * v_m * m_k + g.CD0 * V_l * v_l);

        out.F += F;
        out.M += (r_k - r_cm).cross(F);
    }

    return out;
}

Wrench rocket_aero(const RocketParams &p, const RocketControls &u, const RocketAirData &ad,
                   const Vector3d &w, double x_cm, bool burning) {
    Wrench out = rocket_axial_force(p, ad, burning);
    out += rocket_barrowman_normal(p, ad, x_cm);
    out += rocket_crossflow(p, ad, w, x_cm);
    out += rocket_grid_fins(p, u, ad, w, x_cm);
    return out;
}

// --- gravity ---
Vector3d rocket_gravity(double h) {
    const double r =
        EARTH_RADIUS_MEAN / (EARTH_RADIUS_MEAN + std::max(h, -0.5 * EARTH_RADIUS_MEAN));
    return Vector3d(0.0, 0.0, ISA_G0 * r * r);
}

// --- everything ---
RocketTerms rocket_terms(const RocketParams &p, const RocketEnvironment &env,
                         const RocketControls &u, const RocketAir &air, const RocketState &x,
                         double t) {
    RocketTerms o;
    const double m = x.mass();
    const Vector3d w = x.rates();

    o.mp = rocket_mass_properties(p, m);
    o.engine = rocket_engine(p, u, t, m, o.mp.x_cm, air.p);
    o.air = rocket_air_data(x, air);

    if (env.aero) {
        o.aero = rocket_aero(p, u, o.air, w, o.mp.x_cm, o.engine.T > 0.0);
    }

    // variable mass terms (all zero when the engine is off)
    const double m_dot_e = o.engine.m_dot_e;
    const Vector3d &rho_e = o.engine.rho_e;

    if (env.coriolis) {
        o.F_cor = -2.0 * m_dot_e * w.cross(rho_e);
    }

    if (env.jet_damping) {
        o.M_jet = -m_dot_e * rho_e.cross(w.cross(rho_e));
    }

    if (env.J_dot) {
        o.J_dot = rocket_J_dot(p, m, -m_dot_e);
        o.M_J_dot = -o.J_dot * w;
    }

    if (env.gravity) {
        o.g_i = rocket_gravity(-x.pos().z());
    }

    return o;
}

} // namespace BSIM
