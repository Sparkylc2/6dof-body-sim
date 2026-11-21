// the c api (bsim_c.h) and the packing it relies on (bridge.hpp)
#include "bsim_c.h"
#include "atmosphere.hpp"
#include "bridge.hpp"
#include "craft.hpp"
#include "dynamics.hpp"
#include "tables.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>

namespace BSIM::Bridge {
using Eigen::Map;
using Eigen::Matrix3d;
using Eigen::Vector3d;
using Eigen::Vector4d;

// --- fixed wing ---
#define BSIM_COUNT_FIELD(f) +1
constexpr int FIXEDWING_PARAM_COUNT = 0 BSIM_FIXEDWING_PARAM_FIELDS(BSIM_COUNT_FIELD);
#undef BSIM_COUNT_FIELD

// catches a field being added to FixedWingParams but not to the list (or the header's size)
static_assert(sizeof(FixedWingParams) == FIXEDWING_PARAM_COUNT * sizeof(double),
              "FixedWingParams changed, update BSIM_FIXEDWING_PARAM_FIELDS in bridge.hpp");
static_assert(FIXEDWING_HEAD + FIXEDWING_PARAM_COUNT == BSIM_FIXEDWING_NP,
              "BSIM_FIXEDWING_NP in bsim_c.h is out of date");

FixedWingDef unpack_fixedwing(const double *P) {
    FixedWingDef d;
    d.m = P[0];
    d.J = Map<const Vector4d>(P + 1);
    int i = FIXEDWING_HEAD;
#define BSIM_UNPACK_FIELD(f) d.params.f = P[i++];
    BSIM_FIXEDWING_PARAM_FIELDS(BSIM_UNPACK_FIELD)
#undef BSIM_UNPACK_FIELD
    return d;
}

void pack_fixedwing(const FixedWingDef &d, double *P) {
    P[0] = d.m;
    Map<Vector4d>(P + 1) = d.J;
    int i = FIXEDWING_HEAD;
#define BSIM_PACK_FIELD(f) P[i++] = d.params.f;
    BSIM_FIXEDWING_PARAM_FIELDS(BSIM_PACK_FIELD)
#undef BSIM_PACK_FIELD
}

static FixedWingDef def_of(const FixedWing &plane) { return {plane.m, plane.J, plane.params}; }
FixedWingDef zagi() { return def_of(generate_zagi_flying_wing()); }
FixedWingDef sport() { return def_of(generate_sport_drone()); }

// --- spacecraft ---
int spacecraft_np(int n_thrusters) { return BSIM_SPACECRAFT_NP(n_thrusters); }

bool unpack_spacecraft(const double *P, int nP, SpacecraftDef &out) {
    if (nP < BSIM_SPACECRAFT_NP_BASE) {
        return false;
    }
    const double n_raw = P[BSIM_SPACECRAFT_NP_BASE - 1];
    const int n = static_cast<int>(n_raw);
    if (n < 0 || n != n_raw || nP != spacecraft_np(n)) {
        return false;
    }

    SpacecraftParams &p = out.params;
    p.m_empty = P[0];
    p.J_empty = Map<const Matrix3d>(P + 1);
    p.r_tank = Map<const Vector3d>(P + 10);
    p.R_tank = P[13];
    p.wheels.A = Map<const Matrix3d>(P + 14);
    p.wheels.J_w = P[23];
    p.wheels.tau_max = P[24];
    p.wheels.h_max = P[25];
    out.env.gravity = P[26] != 0.0;
    out.env.j2 = P[27] != 0.0;
    out.env.gravity_gradient = P[28] != 0.0;

    p.thrusters.resize(n);
    for (int k = 0; k < n; k++) {
        const double *t = P + BSIM_SPACECRAFT_NP_BASE + BSIM_THRUSTER_NP * k;
        p.thrusters[k].pos = Map<const Vector3d>(t);
        p.thrusters[k].dir = Map<const Vector3d>(t + 3);
        p.thrusters[k].T_max = t[6];
        p.thrusters[k].Isp = t[7];
    }
    return true;
}

void pack_spacecraft(const SpacecraftDef &d, double *P) {
    const SpacecraftParams &p = d.params;
    const int n = static_cast<int>(p.thrusters.size());
    P[0] = p.m_empty;
    Map<Matrix3d>(P + 1) = p.J_empty;
    Map<Vector3d>(P + 10) = p.r_tank;
    P[13] = p.R_tank;
    Map<Matrix3d>(P + 14) = p.wheels.A;
    P[23] = p.wheels.J_w;
    P[24] = p.wheels.tau_max;
    P[25] = p.wheels.h_max;
    P[26] = d.env.gravity ? 1.0 : 0.0;
    P[27] = d.env.j2 ? 1.0 : 0.0;
    P[28] = d.env.gravity_gradient ? 1.0 : 0.0;
    P[29] = n;
    for (int k = 0; k < n; k++) {
        double *t = P + BSIM_SPACECRAFT_NP_BASE + BSIM_THRUSTER_NP * k;
        Map<Vector3d>{t} = p.thrusters[k].pos;
        Map<Vector3d>(t + 3) = p.thrusters[k].dir;
        t[6] = p.thrusters[k].T_max;
        t[7] = p.thrusters[k].Isp;
    }
}

// --- rocket ---
namespace {

constexpr int RK_ENGINE = 9, RK_FINS = 22, RK_ENV = 29, RK_COUNTS = 34;
constexpr int RK_PROPELLANT_N = 6, RK_COMPONENT_N = 2, RK_SLICE_N = 3;
static_assert(RK_COUNTS + 7 == BSIM_ROCKET_NP_HEAD, "BSIM_ROCKET_NP_HEAD in bsim_c.h is wrong");

// the tables in the order they're packed
const Table1D *rocket_tables(const RocketParams &p, int k) {
    const Table1D *t[] = {&p.engine.thrust_curve, &p.drag.friction, &p.drag.pressure,
                          &p.drag.parasitic};
    return t[k];
}
Table1D *rocket_tables(RocketParams &p, int k) {
    return const_cast<Table1D *>(rocket_tables(std::as_const(p), k));
}

int rocket_np(const int counts[7]) {
    return BSIM_ROCKET_NP_HEAD + RK_PROPELLANT_N * counts[0] + RK_COMPONENT_N * counts[1] +
           RK_SLICE_N * counts[2] + 2 * (counts[3] + counts[4] + counts[5] + counts[6]);
}

void rocket_counts(const RocketParams &p, int counts[7]) {
    counts[0] = static_cast<int>(p.propellant.size());
    counts[1] = static_cast<int>(p.components.size());
    counts[2] = static_cast<int>(p.slices.size());
    for (int k = 0; k < 4; k++) {
        counts[3 + k] = static_cast<int>(rocket_tables(p, k)->x.size());
    }
}

// a whole, non negative number stored as a double
int whole(double v, const char *what) {
    if (!(v >= 0.0) || v != std::floor(v) || v > 1e6) {
        throw std::invalid_argument(std::string(what) + " has to be a whole number >= 0");
    }
    return static_cast<int>(v);
}

} // namespace

int rocket_np(const RocketDef &d) {
    int counts[7];
    rocket_counts(d.params, counts);
    return rocket_np(counts);
}

std::vector<double> pack_rocket(const RocketDef &d) {
    const RocketParams &p = d.params;
    int counts[7];
    rocket_counts(p, counts);
    std::vector<double> P(rocket_np(counts));

    P[0] = p.m_dry;
    P[1] = p.x_dry;
    P[2] = p.Ix_dry;
    P[3] = p.It_dry;
    P[4] = p.d;
    P[5] = p.A_base;
    P[6] = p.A_motor;
    P[7] = p.C_dc;
    P[8] = p.CD_flat;

    const RocketEngine &e = p.engine;
    double *E = P.data() + RK_ENGINE;
    E[0] = e.kind == RocketEngine::Kind::ThrustCurve ? 1.0 : 0.0;
    E[1] = e.c_eff;
    E[2] = e.p_test;
    E[3] = e.T_vac_max;
    E[4] = e.Isp_vac;
    E[5] = e.throttle_min;
    E[6] = e.A_e;
    E[7] = e.x_G;
    E[8] = e.L_n;
    E[9] = e.m_n;
    E[10] = e.l_n;
    E[11] = e.I_na;
    E[12] = e.I_nt;

    double *F = P.data() + RK_FINS;
    F[0] = p.fins.n;
    F[1] = p.fins.x;
    F[2] = p.fins.r;
    F[3] = p.fins.phi0;
    F[4] = p.fins.S;
    F[5] = p.fins.CN_alpha;
    F[6] = p.fins.CD0;

    double *V = P.data() + RK_ENV;
    V[0] = d.env.gravity;
    V[1] = d.env.aero;
    V[2] = d.env.jet_damping;
    V[3] = d.env.J_dot;
    V[4] = d.env.coriolis;

    std::copy(counts, counts + 7, P.begin() + RK_COUNTS);

    double *o = P.data() + BSIM_ROCKET_NP_HEAD;
    for (const PropellantStore &s : p.propellant) {
        *o++ = s.kind == PropellantStore::Kind::Solid ? 1.0 : 0.0;
        *o++ = s.share;
        *o++ = s.x;
        *o++ = s.R;
        *o++ = s.rho;
        *o++ = s.L;
    }
    for (const BodyComponent &c : p.components) {
        *o++ = c.CN_alpha;
        *o++ = c.x_cp;
    }
    for (const BodySlice &s : p.slices) {
        *o++ = s.x;
        *o++ = s.d;
        *o++ = s.dx;
    }
    for (int k = 0; k < 4; k++) {
        const Table1D &t = *rocket_tables(p, k);
        o = std::copy(t.x.begin(), t.x.end(), o);
        o = std::copy(t.y.begin(), t.y.end(), o);
    }
    return P;
}

RocketDef unpack_rocket(const double *P, int nP) {
    if (nP < BSIM_ROCKET_NP_HEAD) {
        throw std::invalid_argument("a packed rocket has at least " +
                                    std::to_string(BSIM_ROCKET_NP_HEAD) + " entries");
    }
    const char *count_names[7] = {"n_propellant", "n_components", "n_slices",   "n_thrust_curve",
                                  "n_friction",   "n_pressure",   "n_parasitic"};
    int counts[7];
    for (int k = 0; k < 7; k++) {
        counts[k] = whole(P[RK_COUNTS + k], count_names[k]);
    }
    if (nP != rocket_np(counts)) {
        throw std::invalid_argument(
            "P has " + std::to_string(nP) + " entries but its counts say it should have " +
            std::to_string(rocket_np(counts)) + " (is it a packed rocket?)");
    }

    RocketDef d;
    RocketParams &p = d.params;
    p.m_dry = P[0];
    p.x_dry = P[1];
    p.Ix_dry = P[2];
    p.It_dry = P[3];
    p.d = P[4];
    p.A_base = P[5];
    p.A_motor = P[6];
    p.C_dc = P[7];
    p.CD_flat = P[8];

    const double *E = P + RK_ENGINE;
    if (E[0] != 0.0 && E[0] != 1.0) {
        throw std::invalid_argument("engine kind has to be 0 (throttleable) or 1 (thrust curve)");
    }
    RocketEngine &e = p.engine;
    e.kind = E[0] == 1.0 ? RocketEngine::Kind::ThrustCurve : RocketEngine::Kind::Throttleable;
    e.c_eff = E[1];
    e.p_test = E[2];
    e.T_vac_max = E[3];
    e.Isp_vac = E[4];
    e.throttle_min = E[5];
    e.A_e = E[6];
    e.x_G = E[7];
    e.L_n = E[8];
    e.m_n = E[9];
    e.l_n = E[10];
    e.I_na = E[11];
    e.I_nt = E[12];

    const double *F = P + RK_FINS;
    p.fins.n = whole(F[0], "fins.n");
    p.fins.x = F[1];
    p.fins.r = F[2];
    p.fins.phi0 = F[3];
    p.fins.S = F[4];
    p.fins.CN_alpha = F[5];
    p.fins.CD0 = F[6];

    const double *V = P + RK_ENV;
    d.env.gravity = V[0] != 0.0;
    d.env.aero = V[1] != 0.0;
    d.env.jet_damping = V[2] != 0.0;
    d.env.J_dot = V[3] != 0.0;
    d.env.coriolis = V[4] != 0.0;

    const double *in = P + BSIM_ROCKET_NP_HEAD;
    p.propellant.resize(counts[0]);
    for (PropellantStore &s : p.propellant) {
        if (in[0] != 0.0 && in[0] != 1.0) {
            throw std::invalid_argument("propellant kind has to be 0 (liquid) or 1 (solid)");
        }
        s.kind = in[0] == 1.0 ? PropellantStore::Kind::Solid : PropellantStore::Kind::Liquid;
        s.share = in[1];
        s.x = in[2];
        s.R = in[3];
        s.rho = in[4];
        s.L = in[5];
        in += RK_PROPELLANT_N;
    }
    p.components.resize(counts[1]);
    for (BodyComponent &c : p.components) {
        c.CN_alpha = in[0];
        c.x_cp = in[1];
        in += RK_COMPONENT_N;
    }
    p.slices.resize(counts[2]);
    for (BodySlice &s : p.slices) {
        s.x = in[0];
        s.d = in[1];
        s.dx = in[2];
        in += RK_SLICE_N;
    }
    const char *table_names[4] = {"engine.thrust_curve", "drag.friction", "drag.pressure",
                                  "drag.parasitic"};
    for (int k = 0; k < 4; k++) {
        const int n = counts[3 + k];
        if (n == 0) {
            continue;
        }
        try {
            *rocket_tables(p, k) = Table1D({in, in + n}, {in + n, in + 2 * n});
        } catch (const std::invalid_argument &err) {
            throw std::invalid_argument(std::string(table_names[k]) + ": " + err.what());
        }
        in += 2 * n;
    }
    return d;
}

RocketExample rocket_example() {
    const Rocket r = generate_example_rocket();
    RocketExample ex;
    ex.def = {r.params, r.env};
    ex.x0 = r.state;
    ex.fin = r.fins.empty() ? Actuator{} : r.fins.front();
    ex.tvc = r.tvc;
    return ex;
}

RocketControls unpack_rocket_controls(const double *u, int n_fins, bool with_motion) {
    RocketControls c;
    c.gimbal_p = u[0];
    c.gimbal_y = u[1];
    c.throttle = u[2];
    c.t_ignition = u[3];
    c.grid_fins_deployed = u[4] != 0.0;
    c.fins.assign(u + 5, u + 5 + n_fins);
    if (with_motion) {
        const double *g = u + 5 + n_fins;
        c.gimbal_p_dot = g[0];
        c.gimbal_y_dot = g[1];
        c.gimbal_p_ddot = g[2];
        c.gimbal_y_ddot = g[3];
    }
    return c;
}

RocketAir unpack_rocket_air(const double *air) {
    RocketAir a;
    a.rho = air[0];
    a.a = air[1];
    a.p = air[2];
    a.wind_i = Map<const Vector3d>(air + 3);
    a.gust_b = Map<const Vector3d>(air + 6);
    return a;
}

// --- tvc ---
#define BSIM_COUNT_FIELD(f) +1
constexpr int TVC_PARAM_COUNT = 0 BSIM_TVC_FIELDS(BSIM_COUNT_FIELD);
#undef BSIM_COUNT_FIELD
static_assert(sizeof(TvcParams) == TVC_PARAM_COUNT * sizeof(double),
              "TvcParams changed, update BSIM_TVC_FIELDS in bridge.hpp");
static_assert(TVC_PARAM_COUNT == BSIM_TVC_NP, "BSIM_TVC_NP in bsim_c.h is out of date");
static_assert(TVC_NX == BSIM_TVC_NX, "BSIM_TVC_NX in bsim_c.h is out of date");

TvcParams unpack_tvc(const double *P) {
    TvcParams p;
    int i = 0;
#define BSIM_UNPACK_FIELD(f) p.f = P[i++];
    BSIM_TVC_FIELDS(BSIM_UNPACK_FIELD)
#undef BSIM_UNPACK_FIELD
    return p;
}

void pack_tvc(const TvcParams &p, double *P) {
    int i = 0;
#define BSIM_PACK_FIELD(f) P[i++] = p.f;
    BSIM_TVC_FIELDS(BSIM_PACK_FIELD)
#undef BSIM_PACK_FIELD
}

// --- air ---
AirState unpack_air(const double *air) {
    AirState a;
    a.rho = air[0];
    a.a = air[1];
    a.wind_i = Map<const Vector3d>(air + 2);
    a.gust_b = Map<const Vector3d>(air + 5);
    return a;
}

} // namespace BSIM::Bridge

// ---- the c api ----
using namespace BSIM;
using namespace BSIM::Bridge;
using Eigen::Map;
using Eigen::Matrix;
using Eigen::Vector3d;
using Eigen::Vector4d;

namespace {

AircraftState aircraft_state(const double *x) {
    AircraftState s;
    s.state = Map<const StateVector<13>>(x);
    return s;
}

FixedWingControls fixedwing_controls(const double *u) { return {u[0], u[1], u[2], u[3]}; }

// x_dot for the fixed wing, with the wrench it came from
StateVector<13> fixedwing_xdot(const FixedWingDef &d, const AircraftState &s,
                               const FixedWingControls &u, const AirState &air, Wrench &w) {
    w = fixed_wing_wrench(d.m, d.params, s, air, u);
    return aircraft_derivatives(d.m, d.J, s, w.F, w.M);
}

} // namespace

extern "C" {

void bsim_aircraft_xdot(double m, const double J[4], const double x[13], const double F[3],
                        const double M[3], double xdot[13]) {
    Map<StateVector<13>>{xdot} =
        aircraft_derivatives(m, Map<const Vector4d>(J), aircraft_state(x), Map<const Vector3d>(F),
                             Map<const Vector3d>(M));
}

void bsim_fixedwing_wrench(const double P[47], const double x[13], const double u[4],
                           const double air[8], double F[3], double M[3]) {
    const FixedWingDef d = unpack_fixedwing(P);
    const Wrench w =
        fixed_wing_wrench(d.m, d.params, aircraft_state(x), unpack_air(air), fixedwing_controls(u));
    Map<Vector3d>{F} = w.F;
    Map<Vector3d>{M} = w.M;
}

void bsim_fixedwing_xdot(const double P[47], const double x[13], const double u[4],
                         const double air[8], double xdot[13], double F[3], double M[3]) {
    Wrench w;
    Map<StateVector<13>>{xdot} = fixedwing_xdot(unpack_fixedwing(P), aircraft_state(x),
                                                fixedwing_controls(u), unpack_air(air), w);
    Map<Vector3d>{F} = w.F;
    Map<Vector3d>{M} = w.M;
}

void bsim_fixedwing_jac(const double P[47], const double x[13], const double u[4],
                        const double air[8], double A[169], double B[52]) {
    const FixedWingDef d = unpack_fixedwing(P);
    const AirState a = unpack_air(air);
    const StateVector<13> x0 = Map<const StateVector<13>>(x);
    const Eigen::Vector4d u0 = Map<const Vector4d>(u);
    Wrench w;

    auto f = [&](const StateVector<13> &xs, const Vector4d &us) {
        AircraftState s;
        s.state = xs;
        return fixedwing_xdot(d, s, fixedwing_controls(us.data()), a, w);
    };
    // step scaled to the size of each entry
    auto step = [](double v) { return 1e-6 * std::max(1.0, std::abs(v)); };

    Map<Matrix<double, 13, 13>> Am(A);
    for (int i = 0; i < 13; i++) {
        const double h = step(x0(i));
        StateVector<13> xp = x0, xm = x0;
        xp(i) += h;
        xm(i) -= h;
        Am.col(i) = (f(xp, u0) - f(xm, u0)) / (2.0 * h);
    }

    Map<Matrix<double, 13, 4>> Bm(B);
    for (int i = 0; i < 4; i++) {
        const double h = step(u0(i));
        Vector4d up = u0, um = u0;
        up(i) += h;
        um(i) -= h;
        Bm.col(i) = (f(x0, up) - f(x0, um)) / (2.0 * h);
    }
}

void bsim_fixedwing_zagi(double P[47]) { pack_fixedwing(zagi(), P); }
void bsim_fixedwing_sport(double P[47]) { pack_fixedwing(sport(), P); }

int bsim_spacecraft_xdot(const double *P, int nP, const double x[17], const double *throttle,
                         const double wheel_torque[3], double xdot[17]) {
    SpacecraftDef d;
    if (!unpack_spacecraft(P, nP, d)) {
        return -1;
    }
    SpacecraftControls u;
    u.throttle.assign(throttle, throttle + d.params.thrusters.size());
    u.wheel_torque = Map<const Vector3d>(wheel_torque);

    SpacecraftState s;
    s.state = Map<const StateVector<17>>(x);
    Map<StateVector<17>>{xdot} = spacecraft_derivatives(d.params, d.env, u, s);
    return 0;
}

int bsim_spacecraft_cubesat(double *P, int nP_max, double x0[17]) {
    const Spacecraft s = generate_cubesat();
    const int nP = spacecraft_np(static_cast<int>(s.params.thrusters.size()));
    if (nP_max >= nP) {
        pack_spacecraft({s.params, s.env}, P);
        Map<StateVector<17>>{x0} = s.state.state;
    }
    return nP;
}

int bsim_rocket_xdot(const double *P, int nP, const double x[14], const double *u, int nU,
                     const double air[9], double t, double xdot[14], double terms[41]) {
    RocketDef d;
    try {
        d = unpack_rocket(P, nP);
    } catch (const std::exception &) {
        std::fill(xdot, xdot + BSIM_ROCKET_NX, std::nan(""));
        std::fill(terms, terms + BSIM_ROCKET_NTERMS, std::nan(""));
        return -1;
    }
    const int n_fins = d.params.fins.n;
    if (nU != BSIM_ROCKET_NU(n_fins) && nU != BSIM_ROCKET_NU_STILL(n_fins)) {
        std::fill(xdot, xdot + BSIM_ROCKET_NX, std::nan(""));
        std::fill(terms, terms + BSIM_ROCKET_NTERMS, std::nan(""));
        return -1;
    }

    RocketState s;
    s.state = Map<const StateVector<14>>(x);
    RocketTerms T;
    Map<StateVector<14>>{xdot} = rocket_derivatives(
        d.params, d.env, unpack_rocket_controls(u, n_fins, nU == BSIM_ROCKET_NU(n_fins)),
        unpack_rocket_air(air), s, t, &T);

    terms[0] = T.mp.m;
    terms[1] = T.mp.x_cm;
    terms[2] = T.mp.Ix;
    terms[3] = T.mp.It;
    terms[4] = T.engine.T;
    terms[5] = T.engine.m_dot_e;
    terms[6] = T.air.V;
    terms[7] = T.air.mach;
    terms[8] = T.air.qbar;
    terms[9] = T.air.alpha_t;
    Map<Vector3d>(terms + 10) = T.F_body();
    Map<Vector3d>(terms + 13) = T.M_body();
    Map<Vector3d>(terms + 16) = T.aero.F;
    Map<Vector3d>(terms + 19) = T.aero.M;
    Map<Vector3d>(terms + 22) = T.g_i;
    Map<Vector3d>(terms + 25) = T.engine.n_T;
    Map<Vector3d>(terms + 28) = T.engine.rho_e;
    Map<Vector3d>(terms + 31) = T.air.v_r;
    Map<Vector3d>(terms + 34) = T.gimbal.M();
    Map<Eigen::Vector2d>(terms + 37) = T.hinge;
    Map<Eigen::Vector2d>(terms + 39) = T.nozzle_load;
    return 0;
}

int bsim_rocket_example(double *P, int nP_max, double x0[14]) {
    const RocketExample ex = rocket_example();
    const std::vector<double> Pv = pack_rocket(ex.def);
    const int nP = static_cast<int>(Pv.size());
    if (nP_max >= nP) {
        std::copy(Pv.begin(), Pv.end(), P);
        Map<StateVector<14>>{x0} = ex.x0.state;
    }
    return nP;
}

void bsim_isa_rocket_air(double h, const double wind_i[3], const double gust_b[3], double air[9]) {
    const Atmosphere atm = isa(h);
    air[0] = atm.rho;
    air[1] = atm.a;
    air[2] = atm.p;
    std::copy(wind_i, wind_i + 3, air + 3);
    std::copy(gust_b, gust_b + 3, air + 6);
}

} // extern "C"

// (c++ helpers for the tvc functions below, outside the extern "C" block)
namespace {

TvcInputs tvc_inputs(const double *delta_cmd, const double *load, double thrust) {
    TvcInputs in;
    in.delta_cmd = Map<const Eigen::Vector2d>(delta_cmd);
    in.load = Map<const Eigen::Vector2d>(load);
    in.thrust = thrust;
    return in;
}

void tvc_outputs(const TvcState &s, const TvcOutputs &o, double *out) {
    Map<Eigen::Vector4d>{out} = s.segment<4>(0); // delta, delta_dot
    Map<Eigen::Vector2d>(out + 4) = o.delta_ddot;
    Map<Eigen::Vector4d>(out + 6) = o.l;
    Map<Eigen::Vector4d>(out + 10) = o.l_cmd;
    Map<Eigen::Vector4d>(out + 14) = o.F;
}

} // namespace

extern "C" {

void bsim_tvc_ember(double P[23]) { pack_tvc(rocket_example().tvc, P); }

void bsim_tvc_step(const double P[23], const double s0[12], const double delta_cmd[2],
                   const double load[2], double thrust, double dt, double s1[12], double out[18]) {
    const TvcState s = Map<const TvcState>(s0);
    TvcOutputs o;
    const TvcState next = tvc_step(unpack_tvc(P), s, tvc_inputs(delta_cmd, load, thrust), dt, &o);
    Map<TvcState>{s1} = next;
    tvc_outputs(next, o, out);
}

void bsim_tvc_xdot(const double P[23], const double s[12], const double delta_cmd[2],
                   const double load[2], double thrust, double sdot[12], double out[18]) {
    const TvcState x = Map<const TvcState>(s);
    TvcOutputs o;
    Map<TvcState>{sdot} =
        tvc_derivatives(unpack_tvc(P), x, tvc_inputs(delta_cmd, load, thrust), &o);
    tvc_outputs(x, o, out);
}

void bsim_tvc_kinematics(const double P[23], const double delta[2], double l[4], double J[8]) {
    const TvcKinematics k = tvc_kinematics(unpack_tvc(P), Map<const Eigen::Vector2d>(delta));
    Map<Eigen::Vector4d>{l} = k.l;
    Map<Matrix42d>{J} = k.J;
}

double bsim_tvc_forward_kinematics(const double P[23], const double l[4], const double delta0[2],
                                   double delta[2]) {
    double res = 0.0;
    Map<Eigen::Vector2d>{delta} = tvc_forward_kinematics(
        unpack_tvc(P), Map<const Eigen::Vector4d>(l), Map<const Eigen::Vector2d>(delta0), 5, &res);
    return res;
}

void bsim_isa(double h, double out[4]) {
    const Atmosphere atm = isa(h);
    out[0] = atm.T;
    out[1] = atm.p;
    out[2] = atm.rho;
    out[3] = atm.a;
}

void bsim_isa_air(double h, const double wind_i[3], const double gust_b[3], double air[8]) {
    const Atmosphere atm = isa(h);
    air[0] = atm.rho;
    air[1] = atm.a;
    std::copy(wind_i, wind_i + 3, air + 2);
    std::copy(gust_b, gust_b + 3, air + 5);
}

double bsim_lookup1(const double *x, const double *y, int n, double xq) {
    const auto N = static_cast<std::size_t>(n);
    return lookup1({x, N}, {y, N}, xq);
}

double bsim_lookup2(const double *x, int nx, const double *y, int ny, const double *z, double xq,
                    double yq) {
    const auto NX = static_cast<std::size_t>(nx), NY = static_cast<std::size_t>(ny);
    return lookup2({x, NX}, {y, NY}, {z, NX * NY}, xq, yq);
}

} // extern "C"
