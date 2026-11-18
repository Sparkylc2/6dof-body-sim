// bsim_mex: one mex for everything, the first argument says what to do
// (the +bsim package wraps these, so normally nothing calls this directly)
//
//   P          = bsim_mex('fixedwing_params', 'zagi' | 'sport')        struct
//   P          = bsim_mex('fixedwing_pack', S)                          struct -> packed
//   S          = bsim_mex('fixedwing_unpack', P)                        packed -> struct
//   [F, M]     = bsim_mex('fixedwing_wrench', P, x, u, air)
//   [xd, F, M] = bsim_mex('fixedwing_xdot', P, x, u, air)
//   [A, B]     = bsim_mex('fixedwing_jac', P, x, u, air)
//   xd         = bsim_mex('aircraft_xdot', m, J, x, F, M)
//   [S, x0]    = bsim_mex('spacecraft_params', 'cubesat')
//   P          = bsim_mex('spacecraft_pack', S)
//   S          = bsim_mex('spacecraft_unpack', P)
//   xd         = bsim_mex('spacecraft_xdot', P, x, throttle, wheel_torque)
//   [S, x0, act] = bsim_mex('rocket_params', 'example')                 (act: actuator numbers)
//   P          = bsim_mex('rocket_pack', S)
//   S          = bsim_mex('rocket_unpack', P)
//   [xd, terms, tv] = bsim_mex('rocket_xdot', P, x, u, air, t)          (terms struct, tv flat)
//   air        = bsim_mex('isa_rocket_air', h, wind_i, gust_b)
//   [T, p, rho, a] = bsim_mex('isa', h)                                 (h can be any size)
//   air        = bsim_mex('isa_air', h, wind_i, gust_b)
//   yq         = bsim_mex('lookup1', x, y, xq)                          (xq can be any size)
//   zq         = bsim_mex('lookup2', x, y, Z, xq, yq)                   (Z is numel(x) by numel(y))
//
// wherever P is taken it can be the struct or the packed vector (packed skips the conversion, so
// it's the one to use inside loops)
#include "bridge.hpp"
#include "bsim_c.h"
#include "mex.h"
#include "tables.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <vector>

using namespace BSIM;
using namespace BSIM::Bridge;

namespace {

// errors are thrown as exceptions and only turned into a matlab error at the very top, after
// everything on the c++ side has been destroyed (mexErrMsgIdAndTxt doesn't unwind the stack)
struct Error : std::runtime_error {
    using std::runtime_error::runtime_error;
};

// --- reading arguments ---
std::string get_string(const mxArray *a, const char *what) {
    if (!mxIsChar(a)) {
        throw Error(std::string(what) + " has to be a string");
    }
    char *s = mxArrayToString(a);
    std::string out(s);
    mxFree(s);
    return out;
}

// a real double array with exactly n entries (any shape)
const double *get_doubles(const mxArray *a, std::size_t n, const char *what) {
    if (!mxIsDouble(a) || mxIsComplex(a) || mxIsSparse(a)) {
        throw Error(std::string(what) + " has to be a real double array");
    }
    if (n != 0 && mxGetNumberOfElements(a) != n) {
        throw Error(std::string(what) + " has to have " + std::to_string(n) + " entries, not " +
                    std::to_string(mxGetNumberOfElements(a)));
    }
    return mxGetDoubles(a);
}

double get_scalar(const mxArray *a, const char *what) {
    if (mxIsLogical(a) && mxGetNumberOfElements(a) == 1) {
        return mxIsLogicalScalarTrue(a) ? 1.0 : 0.0;
    }
    return *get_doubles(a, 1, what);
}

const mxArray *get_field(const mxArray *s, const char *name, mwIndex i = 0) {
    const mxArray *f = mxGetField(s, i, name);
    if (!f) {
        throw Error(std::string("missing field '") + name + "'");
    }
    return f;
}

double field_scalar(const mxArray *s, const char *name, mwIndex i = 0) {
    return get_scalar(get_field(s, name, i), name);
}

const double *field_doubles(const mxArray *s, const char *name, std::size_t n, mwIndex i = 0) {
    return get_doubles(get_field(s, name, i), n, name);
}

void check_args(int nrhs, int needed, const std::string &cmd) {
    if (nrhs != needed) {
        throw Error(cmd + " takes " + std::to_string(needed - 1) + " arguments");
    }
}

// --- making outputs ---
mxArray *column(const double *v, std::size_t n) {
    mxArray *a = mxCreateDoubleMatrix(n, 1, mxREAL);
    std::memcpy(mxGetDoubles(a), v, n * sizeof(double));
    return a;
}

mxArray *matrix(const double *v, std::size_t rows, std::size_t cols) {
    mxArray *a = mxCreateDoubleMatrix(rows, cols, mxREAL);
    std::memcpy(mxGetDoubles(a), v, rows * cols * sizeof(double));
    return a;
}

mxArray *scalar(double v) { return mxCreateDoubleScalar(v); }

// --- fixed wing structs ---
#define BSIM_FIELD_NAME(f) #f,
const char *FIXEDWING_FIELDS[] = {"m", "J", BSIM_FIXEDWING_PARAM_FIELDS(BSIM_FIELD_NAME)};
#undef BSIM_FIELD_NAME
constexpr int N_FIXEDWING_FIELDS = sizeof(FIXEDWING_FIELDS) / sizeof(FIXEDWING_FIELDS[0]);

mxArray *fixedwing_struct(const FixedWingDef &d) {
    mxArray *s = mxCreateStructMatrix(1, 1, N_FIXEDWING_FIELDS, FIXEDWING_FIELDS);
    mxSetField(s, 0, "m", scalar(d.m));
    mxSetField(s, 0, "J", column(d.J.data(), 4));
#define BSIM_SET_FIELD(f) mxSetField(s, 0, #f, scalar(d.params.f));
    BSIM_FIXEDWING_PARAM_FIELDS(BSIM_SET_FIELD)
#undef BSIM_SET_FIELD
    return s;
}

// every field has to be there, and nothing else (so a typo doesn't silently get ignored)
FixedWingDef fixedwing_from_struct(const mxArray *s) {
    const int nf = mxGetNumberOfFields(s);
    for (int k = 0; k < nf; k++) {
        const char *name = mxGetFieldNameByNumber(s, k);
        bool known = false;
        for (const char *f : FIXEDWING_FIELDS) {
            known = known || std::strcmp(name, f) == 0;
        }
        if (!known) {
            throw Error(std::string("unknown fixed wing parameter '") + name + "'");
        }
    }
    FixedWingDef d;
    d.m = field_scalar(s, "m");
    d.J = Eigen::Map<const Eigen::Vector4d>(field_doubles(s, "J", 4));
#define BSIM_GET_FIELD(f) d.params.f = field_scalar(s, #f);
    BSIM_FIXEDWING_PARAM_FIELDS(BSIM_GET_FIELD)
#undef BSIM_GET_FIELD
    return d;
}

// P as either the struct or the packed vector, packed into buf if it was a struct
const double *fixedwing_P(const mxArray *a, double (&buf)[BSIM_FIXEDWING_NP]) {
    if (mxIsStruct(a)) {
        pack_fixedwing(fixedwing_from_struct(a), buf);
        return buf;
    }
    return get_doubles(a, BSIM_FIXEDWING_NP, "P");
}

// --- spacecraft structs ---
mxArray *spacecraft_struct(const SpacecraftDef &d) {
    const SpacecraftParams &p = d.params;

    const char *wheel_fields[] = {"A", "J_w", "tau_max", "h_max"};
    mxArray *wheels = mxCreateStructMatrix(1, 1, 4, wheel_fields);
    mxSetField(wheels, 0, "A", matrix(p.wheels.A.data(), 3, 3));
    mxSetField(wheels, 0, "J_w", scalar(p.wheels.J_w));
    mxSetField(wheels, 0, "tau_max", scalar(p.wheels.tau_max));
    mxSetField(wheels, 0, "h_max", scalar(p.wheels.h_max));

    const char *env_fields[] = {"gravity", "j2", "gravity_gradient"};
    mxArray *env = mxCreateStructMatrix(1, 1, 3, env_fields);
    mxSetField(env, 0, "gravity", mxCreateLogicalScalar(d.env.gravity));
    mxSetField(env, 0, "j2", mxCreateLogicalScalar(d.env.j2));
    mxSetField(env, 0, "gravity_gradient", mxCreateLogicalScalar(d.env.gravity_gradient));

    const char *thr_fields[] = {"pos", "dir", "T_max", "Isp"};
    const auto n = static_cast<mwSize>(p.thrusters.size());
    mxArray *thr = mxCreateStructMatrix(1, n, 4, thr_fields);
    for (mwIndex k = 0; k < n; k++) {
        mxSetField(thr, k, "pos", column(p.thrusters[k].pos.data(), 3));
        mxSetField(thr, k, "dir", column(p.thrusters[k].dir.data(), 3));
        mxSetField(thr, k, "T_max", scalar(p.thrusters[k].T_max));
        mxSetField(thr, k, "Isp", scalar(p.thrusters[k].Isp));
    }

    const char *fields[] = {"m_empty", "J_empty", "r_tank", "R_tank", "wheels", "env", "thrusters"};
    mxArray *s = mxCreateStructMatrix(1, 1, 7, fields);
    mxSetField(s, 0, "m_empty", scalar(p.m_empty));
    mxSetField(s, 0, "J_empty", matrix(p.J_empty.data(), 3, 3));
    mxSetField(s, 0, "r_tank", column(p.r_tank.data(), 3));
    mxSetField(s, 0, "R_tank", scalar(p.R_tank));
    mxSetField(s, 0, "wheels", wheels);
    mxSetField(s, 0, "env", env);
    mxSetField(s, 0, "thrusters", thr);
    return s;
}

SpacecraftDef spacecraft_from_struct(const mxArray *s) {
    using Eigen::Map;
    using Eigen::Matrix3d;
    SpacecraftDef d;
    SpacecraftParams &p = d.params;
    p.m_empty = field_scalar(s, "m_empty");
    p.J_empty = Map<const Matrix3d>(field_doubles(s, "J_empty", 9));
    p.r_tank = Map<const Vector3d>(field_doubles(s, "r_tank", 3));
    p.R_tank = field_scalar(s, "R_tank");

    const mxArray *w = get_field(s, "wheels");
    p.wheels.A = Map<const Matrix3d>(field_doubles(w, "A", 9));
    p.wheels.J_w = field_scalar(w, "J_w");
    p.wheels.tau_max = field_scalar(w, "tau_max");
    p.wheels.h_max = field_scalar(w, "h_max");

    const mxArray *env = get_field(s, "env");
    d.env.gravity = field_scalar(env, "gravity") != 0.0;
    d.env.j2 = field_scalar(env, "j2") != 0.0;
    d.env.gravity_gradient = field_scalar(env, "gravity_gradient") != 0.0;

    const mxArray *thr = get_field(s, "thrusters");
    const std::size_t n = mxIsEmpty(thr) ? 0 : mxGetNumberOfElements(thr);
    if (n > 0 && !mxIsStruct(thr)) {
        throw Error("thrusters has to be a struct array");
    }
    p.thrusters.resize(n);
    for (std::size_t k = 0; k < n; k++) {
        p.thrusters[k].pos = Map<const Vector3d>(field_doubles(thr, "pos", 3, k));
        p.thrusters[k].dir = Map<const Vector3d>(field_doubles(thr, "dir", 3, k));
        p.thrusters[k].T_max = field_scalar(thr, "T_max", k);
        p.thrusters[k].Isp = field_scalar(thr, "Isp", k);
    }
    return d;
}

std::vector<double> pack(const SpacecraftDef &d) {
    std::vector<double> P(spacecraft_np(static_cast<int>(d.params.thrusters.size())));
    pack_spacecraft(d, P.data());
    return P;
}

// P as either the struct or the packed vector
std::vector<double> spacecraft_P(const mxArray *a) {
    if (mxIsStruct(a)) {
        return pack(spacecraft_from_struct(a));
    }
    const double *P = get_doubles(a, 0, "P");
    return {P, P + mxGetNumberOfElements(a)};
}

// --- rocket structs ---
// a typo'd field is an error rather than silently ignored
void check_fields(const mxArray *s, std::initializer_list<const char *> allowed, const char *what) {
    if (!mxIsStruct(s)) {
        throw Error(std::string(what) + " has to be a struct");
    }
    const int nf = mxGetNumberOfFields(s);
    for (int k = 0; k < nf; k++) {
        const char *name = mxGetFieldNameByNumber(s, k);
        bool known = false;
        for (const char *f : allowed) {
            known = known || std::strcmp(name, f) == 0;
        }
        if (!known) {
            throw Error(std::string("unknown field '") + name + "' in " + what);
        }
    }
}

mxArray *struct_of(std::initializer_list<const char *> fields, mwSize n = 1) {
    std::vector<const char *> f(fields);
    return mxCreateStructMatrix(1, n, static_cast<int>(f.size()), f.data());
}

// n rows by cols columns, rows[i * cols + j] row major in, column major out
mxArray *rows_matrix(const std::vector<double> &rows, std::size_t cols) {
    const std::size_t n = rows.size() / cols;
    mxArray *a = mxCreateDoubleMatrix(n, cols, mxREAL);
    double *o = mxGetDoubles(a);
    for (std::size_t i = 0; i < n; i++) {
        for (std::size_t j = 0; j < cols; j++) {
            o[i + n * j] = rows[i * cols + j];
        }
    }
    return a;
}

// an n by cols matrix (or empty), as row major
std::vector<double> matrix_rows(const mxArray *a, std::size_t cols, const char *what) {
    if (mxIsEmpty(a)) {
        return {};
    }
    const double *v = get_doubles(a, 0, what);
    if (mxGetNumberOfDimensions(a) != 2 || mxGetN(a) != cols) {
        throw Error(std::string(what) + " has to have " + std::to_string(cols) + " columns");
    }
    const std::size_t n = mxGetM(a);
    std::vector<double> rows(n * cols);
    for (std::size_t i = 0; i < n; i++) {
        for (std::size_t j = 0; j < cols; j++) {
            rows[i * cols + j] = v[i + n * j];
        }
    }
    return rows;
}

// a table as [x y] columns
mxArray *table_matrix(const Table1D &t) {
    const std::size_t n = t.x.size();
    mxArray *a = mxCreateDoubleMatrix(n, n == 0 ? 0 : 2, mxREAL);
    double *o = mxGetDoubles(a);
    std::copy(t.x.begin(), t.x.end(), o);
    std::copy(t.y.begin(), t.y.end(), o + n);
    return a;
}

Table1D matrix_table(const mxArray *a, const char *what) {
    const std::vector<double> r = matrix_rows(a, 2, what);
    const std::size_t n = r.size() / 2;
    if (n == 0) {
        return {};
    }
    std::vector<double> x(n), y(n);
    for (std::size_t i = 0; i < n; i++) {
        x[i] = r[2 * i];
        y[i] = r[2 * i + 1];
    }
    try {
        return Table1D(std::move(x), std::move(y));
    } catch (const std::invalid_argument &e) {
        throw Error(std::string(what) + ": " + e.what());
    }
}

mxArray *actuator_struct(const Actuator &a) {
    mxArray *s = struct_of({"kind", "tau", "wn", "zeta", "rate_max", "pos_max"});
    const char *kind = a.kind == Actuator::Kind::Ideal        ? "ideal"
                       : a.kind == Actuator::Kind::FirstOrder ? "first_order"
                                                              : "second_order";
    mxSetField(s, 0, "kind", mxCreateString(kind));
    mxSetField(s, 0, "tau", scalar(a.tau));
    mxSetField(s, 0, "wn", scalar(a.wn));
    mxSetField(s, 0, "zeta", scalar(a.zeta));
    mxSetField(s, 0, "rate_max", scalar(a.rate_max));
    mxSetField(s, 0, "pos_max", scalar(a.pos_max));
    return s;
}

mxArray *rocket_struct(const RocketDef &d) {
    const RocketParams &p = d.params;

    const auto n_prop = static_cast<mwSize>(p.propellant.size());
    mxArray *prop = struct_of({"kind", "share", "x", "R", "rho", "L"}, n_prop);
    for (mwIndex k = 0; k < n_prop; k++) {
        const PropellantStore &s = p.propellant[k];
        mxSetField(prop, k, "kind",
                   mxCreateString(s.kind == PropellantStore::Kind::Solid ? "solid" : "liquid"));
        mxSetField(prop, k, "share", scalar(s.share));
        mxSetField(prop, k, "x", scalar(s.x));
        mxSetField(prop, k, "R", scalar(s.R));
        mxSetField(prop, k, "rho", scalar(s.rho));
        mxSetField(prop, k, "L", scalar(s.L));
    }

    const RocketEngine &e = p.engine;
    mxArray *eng = struct_of({"kind", "thrust_curve", "c_eff", "p_test", "T_vac_max", "Isp_vac",
                              "throttle_min", "A_e", "x_G", "L_n"});
    mxSetField(eng, 0, "kind",
               mxCreateString(e.kind == RocketEngine::Kind::ThrustCurve ? "thrust_curve"
                                                                        : "throttleable"));
    mxSetField(eng, 0, "thrust_curve", table_matrix(e.thrust_curve));
    mxSetField(eng, 0, "c_eff", scalar(e.c_eff));
    mxSetField(eng, 0, "p_test", scalar(e.p_test));
    mxSetField(eng, 0, "T_vac_max", scalar(e.T_vac_max));
    mxSetField(eng, 0, "Isp_vac", scalar(e.Isp_vac));
    mxSetField(eng, 0, "throttle_min", scalar(e.throttle_min));
    mxSetField(eng, 0, "A_e", scalar(e.A_e));
    mxSetField(eng, 0, "x_G", scalar(e.x_G));
    mxSetField(eng, 0, "L_n", scalar(e.L_n));

    const auto n_comp = static_cast<mwSize>(p.components.size());
    mxArray *comp = struct_of({"CN_alpha", "x_cp"}, n_comp);
    for (mwIndex k = 0; k < n_comp; k++) {
        mxSetField(comp, k, "CN_alpha", scalar(p.components[k].CN_alpha));
        mxSetField(comp, k, "x_cp", scalar(p.components[k].x_cp));
    }

    std::vector<double> slices;
    for (const BodySlice &s : p.slices) {
        slices.insert(slices.end(), {s.x, s.d, s.dx});
    }

    mxArray *drag = struct_of({"friction", "pressure", "parasitic"});
    mxSetField(drag, 0, "friction", table_matrix(p.drag.friction));
    mxSetField(drag, 0, "pressure", table_matrix(p.drag.pressure));
    mxSetField(drag, 0, "parasitic", table_matrix(p.drag.parasitic));

    mxArray *fins = struct_of({"n", "x", "r", "phi0", "S", "CN_alpha", "CD0"});
    mxSetField(fins, 0, "n", scalar(p.fins.n));
    mxSetField(fins, 0, "x", scalar(p.fins.x));
    mxSetField(fins, 0, "r", scalar(p.fins.r));
    mxSetField(fins, 0, "phi0", scalar(p.fins.phi0));
    mxSetField(fins, 0, "S", scalar(p.fins.S));
    mxSetField(fins, 0, "CN_alpha", scalar(p.fins.CN_alpha));
    mxSetField(fins, 0, "CD0", scalar(p.fins.CD0));

    mxArray *env = struct_of({"gravity", "aero", "jet_damping", "J_dot", "coriolis"});
    mxSetField(env, 0, "gravity", mxCreateLogicalScalar(d.env.gravity));
    mxSetField(env, 0, "aero", mxCreateLogicalScalar(d.env.aero));
    mxSetField(env, 0, "jet_damping", mxCreateLogicalScalar(d.env.jet_damping));
    mxSetField(env, 0, "J_dot", mxCreateLogicalScalar(d.env.J_dot));
    mxSetField(env, 0, "coriolis", mxCreateLogicalScalar(d.env.coriolis));

    mxArray *s =
        struct_of({"m_dry", "x_dry", "Ix_dry", "It_dry", "propellant", "engine", "d", "A_base",
                   "A_motor", "components", "slices", "C_dc", "drag", "CD_flat", "fins", "env"});
    mxSetField(s, 0, "m_dry", scalar(p.m_dry));
    mxSetField(s, 0, "x_dry", scalar(p.x_dry));
    mxSetField(s, 0, "Ix_dry", scalar(p.Ix_dry));
    mxSetField(s, 0, "It_dry", scalar(p.It_dry));
    mxSetField(s, 0, "propellant", prop);
    mxSetField(s, 0, "engine", eng);
    mxSetField(s, 0, "d", scalar(p.d));
    mxSetField(s, 0, "A_base", scalar(p.A_base));
    mxSetField(s, 0, "A_motor", scalar(p.A_motor));
    mxSetField(s, 0, "components", comp);
    mxSetField(s, 0, "slices", rows_matrix(slices, 3));
    mxSetField(s, 0, "C_dc", scalar(p.C_dc));
    mxSetField(s, 0, "drag", drag);
    mxSetField(s, 0, "CD_flat", scalar(p.CD_flat));
    mxSetField(s, 0, "fins", fins);
    mxSetField(s, 0, "env", env);
    return s;
}

// a struct array with n entries (empty is fine), checked for unknown fields
std::size_t struct_array(const mxArray *a, std::initializer_list<const char *> allowed,
                         const char *what) {
    if (mxIsEmpty(a)) {
        return 0;
    }
    check_fields(a, allowed, what);
    return mxGetNumberOfElements(a);
}

RocketDef rocket_from_struct(const mxArray *s) {
    check_fields(s,
                 {"m_dry", "x_dry", "Ix_dry", "It_dry", "propellant", "engine", "d", "A_base",
                  "A_motor", "components", "slices", "C_dc", "drag", "CD_flat", "fins", "env"},
                 "the rocket");
    RocketDef d;
    RocketParams &p = d.params;
    p.m_dry = field_scalar(s, "m_dry");
    p.x_dry = field_scalar(s, "x_dry");
    p.Ix_dry = field_scalar(s, "Ix_dry");
    p.It_dry = field_scalar(s, "It_dry");
    p.d = field_scalar(s, "d");
    p.A_base = field_scalar(s, "A_base");
    p.A_motor = field_scalar(s, "A_motor");
    p.C_dc = field_scalar(s, "C_dc");
    p.CD_flat = field_scalar(s, "CD_flat");

    const mxArray *prop = get_field(s, "propellant");
    const std::size_t n_prop =
        struct_array(prop, {"kind", "share", "x", "R", "rho", "L"}, "propellant");
    p.propellant.resize(n_prop);
    for (std::size_t k = 0; k < n_prop; k++) {
        PropellantStore &st = p.propellant[k];
        const std::string kind = get_string(get_field(prop, "kind", k), "propellant kind");
        if (kind != "liquid" && kind != "solid") {
            throw Error("propellant kind has to be 'liquid' or 'solid', not '" + kind + "'");
        }
        st.kind = kind == "solid" ? PropellantStore::Kind::Solid : PropellantStore::Kind::Liquid;
        st.share = field_scalar(prop, "share", k);
        st.x = field_scalar(prop, "x", k);
        st.R = field_scalar(prop, "R", k);
        st.rho = field_scalar(prop, "rho", k);
        st.L = field_scalar(prop, "L", k);
    }

    const mxArray *eng = get_field(s, "engine");
    check_fields(eng,
                 {"kind", "thrust_curve", "c_eff", "p_test", "T_vac_max", "Isp_vac", "throttle_min",
                  "A_e", "x_G", "L_n"},
                 "engine");
    RocketEngine &e = p.engine;
    const std::string kind = get_string(get_field(eng, "kind"), "engine kind");
    if (kind != "throttleable" && kind != "thrust_curve") {
        throw Error("engine kind has to be 'throttleable' or 'thrust_curve', not '" + kind + "'");
    }
    e.kind =
        kind == "thrust_curve" ? RocketEngine::Kind::ThrustCurve : RocketEngine::Kind::Throttleable;
    e.thrust_curve =
        matrix_table(get_field(eng, "thrust_curve"), "engine.thrust_curve ([t T] columns)");
    e.c_eff = field_scalar(eng, "c_eff");
    e.p_test = field_scalar(eng, "p_test");
    e.T_vac_max = field_scalar(eng, "T_vac_max");
    e.Isp_vac = field_scalar(eng, "Isp_vac");
    e.throttle_min = field_scalar(eng, "throttle_min");
    e.A_e = field_scalar(eng, "A_e");
    e.x_G = field_scalar(eng, "x_G");
    e.L_n = field_scalar(eng, "L_n");

    const mxArray *comp = get_field(s, "components");
    const std::size_t n_comp = struct_array(comp, {"CN_alpha", "x_cp"}, "components");
    p.components.resize(n_comp);
    for (std::size_t k = 0; k < n_comp; k++) {
        p.components[k].CN_alpha = field_scalar(comp, "CN_alpha", k);
        p.components[k].x_cp = field_scalar(comp, "x_cp", k);
    }

    const std::vector<double> slices =
        matrix_rows(get_field(s, "slices"), 3, "slices ([x d dx] columns)");
    p.slices.resize(slices.size() / 3);
    for (std::size_t k = 0; k < p.slices.size(); k++) {
        p.slices[k] = {.x = slices[3 * k], .d = slices[3 * k + 1], .dx = slices[3 * k + 2]};
    }

    const mxArray *drag = get_field(s, "drag");
    check_fields(drag, {"friction", "pressure", "parasitic"}, "drag");
    p.drag.friction = matrix_table(get_field(drag, "friction"), "drag.friction ([mach cd])");
    p.drag.pressure = matrix_table(get_field(drag, "pressure"), "drag.pressure ([mach cd])");
    p.drag.parasitic = matrix_table(get_field(drag, "parasitic"), "drag.parasitic ([mach cd])");

    const mxArray *fins = get_field(s, "fins");
    check_fields(fins, {"n", "x", "r", "phi0", "S", "CN_alpha", "CD0"}, "fins");
    const double n_fins = field_scalar(fins, "n");
    if (!(n_fins >= 0.0) || n_fins != std::floor(n_fins)) {
        throw Error("fins.n has to be a whole number >= 0");
    }
    p.fins.n = static_cast<int>(n_fins);
    p.fins.x = field_scalar(fins, "x");
    p.fins.r = field_scalar(fins, "r");
    p.fins.phi0 = field_scalar(fins, "phi0");
    p.fins.S = field_scalar(fins, "S");
    p.fins.CN_alpha = field_scalar(fins, "CN_alpha");
    p.fins.CD0 = field_scalar(fins, "CD0");

    const mxArray *env = get_field(s, "env");
    check_fields(env, {"gravity", "aero", "jet_damping", "J_dot", "coriolis"}, "env");
    d.env.gravity = field_scalar(env, "gravity") != 0.0;
    d.env.aero = field_scalar(env, "aero") != 0.0;
    d.env.jet_damping = field_scalar(env, "jet_damping") != 0.0;
    d.env.J_dot = field_scalar(env, "J_dot") != 0.0;
    d.env.coriolis = field_scalar(env, "coriolis") != 0.0;
    return d;
}

// P as either the struct or the packed vector (checked either way)
std::vector<double> rocket_P(const mxArray *a) {
    if (mxIsStruct(a)) {
        return pack_rocket(rocket_from_struct(a));
    }
    const double *P = get_doubles(a, 0, "P");
    std::vector<double> v(P, P + mxGetNumberOfElements(a));
    unpack_rocket(v.data(), static_cast<int>(v.size())); // throws if it isn't a rocket
    return v;
}

mxArray *rocket_terms_struct(const double *t) {
    mxArray *s = struct_of({"m", "x_cm", "Ix", "It", "T", "m_dot_e", "V", "mach", "qbar", "alpha_t",
                            "F_body", "M_body", "F_aero", "M_aero", "g_i", "n_T", "rho_e", "v_r"});
    const char *scalars[] = {"m",       "x_cm", "Ix",   "It",   "T",
                             "m_dot_e", "V",    "mach", "qbar", "alpha_t"};
    for (int k = 0; k < 10; k++) {
        mxSetField(s, 0, scalars[k], scalar(t[k]));
    }
    mxSetField(s, 0, "F_body", column(t + 10, 3));
    mxSetField(s, 0, "M_body", column(t + 13, 3));
    mxSetField(s, 0, "F_aero", column(t + 16, 3));
    mxSetField(s, 0, "M_aero", column(t + 19, 3));
    mxSetField(s, 0, "g_i", column(t + 22, 3));
    mxSetField(s, 0, "n_T", column(t + 25, 3));
    mxSetField(s, 0, "rho_e", column(t + 28, 3));
    mxSetField(s, 0, "v_r", column(t + 31, 3));
    return s;
}

// --- commands ---
void dispatch(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 1) {
        throw Error("bsim_mex needs a command, see the top of bsim_mex.cpp");
    }
    const std::string cmd = get_string(prhs[0], "the command");
    (void)nlhs;

    if (cmd == "fixedwing_params") {
        check_args(nrhs, 2, cmd);
        const std::string name = get_string(prhs[1], "the plane");
        if (name == "zagi") {
            plhs[0] = fixedwing_struct(zagi());
        } else if (name == "sport") {
            plhs[0] = fixedwing_struct(sport());
        } else {
            throw Error("no plane called '" + name + "' (there's zagi and sport)");
        }

    } else if (cmd == "fixedwing_pack") {
        check_args(nrhs, 2, cmd);
        double P[BSIM_FIXEDWING_NP];
        plhs[0] = column(fixedwing_P(prhs[1], P), BSIM_FIXEDWING_NP);

    } else if (cmd == "fixedwing_unpack") {
        check_args(nrhs, 2, cmd);
        plhs[0] = fixedwing_struct(unpack_fixedwing(get_doubles(prhs[1], BSIM_FIXEDWING_NP, "P")));

    } else if (cmd == "fixedwing_wrench" || cmd == "fixedwing_xdot" || cmd == "fixedwing_jac") {
        check_args(nrhs, 5, cmd);
        double buf[BSIM_FIXEDWING_NP];
        const double *P = fixedwing_P(prhs[1], buf);
        const double *x = get_doubles(prhs[2], BSIM_AIRCRAFT_NX, "x");
        const double *u = get_doubles(prhs[3], BSIM_FIXEDWING_NU, "u");
        const double *air = get_doubles(prhs[4], BSIM_AIR_N, "air");

        if (cmd == "fixedwing_wrench") {
            double F[3], M[3];
            bsim_fixedwing_wrench(P, x, u, air, F, M);
            plhs[0] = column(F, 3);
            plhs[1] = column(M, 3);
        } else if (cmd == "fixedwing_xdot") {
            double xd[13], F[3], M[3];
            bsim_fixedwing_xdot(P, x, u, air, xd, F, M);
            plhs[0] = column(xd, 13);
            plhs[1] = column(F, 3);
            plhs[2] = column(M, 3);
        } else {
            double A[13 * 13], B[13 * 4];
            bsim_fixedwing_jac(P, x, u, air, A, B);
            plhs[0] = matrix(A, 13, 13);
            plhs[1] = matrix(B, 13, 4);
        }

    } else if (cmd == "aircraft_xdot") {
        check_args(nrhs, 6, cmd);
        double xd[13];
        bsim_aircraft_xdot(get_scalar(prhs[1], "m"), get_doubles(prhs[2], 4, "J"),
                           get_doubles(prhs[3], 13, "x"), get_doubles(prhs[4], 3, "F"),
                           get_doubles(prhs[5], 3, "M"), xd);
        plhs[0] = column(xd, 13);

    } else if (cmd == "spacecraft_params") {
        check_args(nrhs, 2, cmd);
        const std::string name = get_string(prhs[1], "the spacecraft");
        if (name != "cubesat") {
            throw Error("no spacecraft called '" + name + "' (there's cubesat)");
        }
        std::vector<double> P(bsim_spacecraft_cubesat(nullptr, 0, nullptr));
        double x0[17];
        bsim_spacecraft_cubesat(P.data(), static_cast<int>(P.size()), x0);
        SpacecraftDef d;
        unpack_spacecraft(P.data(), static_cast<int>(P.size()), d);
        plhs[0] = spacecraft_struct(d);
        plhs[1] = column(x0, 17);

    } else if (cmd == "spacecraft_pack") {
        check_args(nrhs, 2, cmd);
        const std::vector<double> P = spacecraft_P(prhs[1]);
        plhs[0] = column(P.data(), P.size());

    } else if (cmd == "spacecraft_unpack") {
        check_args(nrhs, 2, cmd);
        const std::vector<double> P = spacecraft_P(prhs[1]);
        SpacecraftDef d;
        if (!unpack_spacecraft(P.data(), static_cast<int>(P.size()), d)) {
            throw Error(
                "P isn't a packed spacecraft (its length doesn't match its thruster count)");
        }
        plhs[0] = spacecraft_struct(d);

    } else if (cmd == "spacecraft_xdot") {
        check_args(nrhs, 5, cmd);
        const std::vector<double> P = spacecraft_P(prhs[1]);
        SpacecraftDef d;
        if (!unpack_spacecraft(P.data(), static_cast<int>(P.size()), d)) {
            throw Error(
                "P isn't a packed spacecraft (its length doesn't match its thruster count)");
        }
        const std::size_t n_thr = d.params.thrusters.size();
        const double *throttle =
            n_thr == 0 ? nullptr : get_doubles(prhs[3], n_thr, "throttle (one per thruster)");
        double xd[17];
        bsim_spacecraft_xdot(P.data(), static_cast<int>(P.size()), get_doubles(prhs[2], 17, "x"),
                             throttle, get_doubles(prhs[4], 3, "wheel_torque"), xd);
        plhs[0] = column(xd, 17);

    } else if (cmd == "rocket_params") {
        check_args(nrhs, 2, cmd);
        const std::string name = get_string(prhs[1], "the rocket");
        if (name != "example") {
            throw Error("no rocket called '" + name + "' (there's example)");
        }
        const RocketExample ex = rocket_example();
        plhs[0] = rocket_struct(ex.def);
        plhs[1] = column(ex.x0.state.data(), BSIM_ROCKET_NX);
        mxArray *act = struct_of({"gimbal", "fin"});
        mxSetField(act, 0, "gimbal", actuator_struct(ex.gimbal));
        mxSetField(act, 0, "fin", actuator_struct(ex.fin));
        plhs[2] = act;

    } else if (cmd == "rocket_pack") {
        check_args(nrhs, 2, cmd);
        const std::vector<double> P = rocket_P(prhs[1]);
        plhs[0] = column(P.data(), P.size());

    } else if (cmd == "rocket_unpack") {
        check_args(nrhs, 2, cmd);
        const std::vector<double> P = rocket_P(prhs[1]);
        plhs[0] = rocket_struct(unpack_rocket(P.data(), static_cast<int>(P.size())));

    } else if (cmd == "rocket_xdot") {
        check_args(nrhs, 6, cmd);
        const std::vector<double> P = rocket_P(prhs[1]);
        const int n_fins = unpack_rocket(P.data(), static_cast<int>(P.size())).params.fins.n;
        const int nU = BSIM_ROCKET_NU(n_fins);
        const std::string u_what = "u ([gimbal_p gimbal_y throttle t_ignition fins_deployed] + " +
                                   std::to_string(n_fins) + " fins)";
        const double *u = get_doubles(prhs[3], nU, u_what.c_str());
        double xd[BSIM_ROCKET_NX], terms[BSIM_ROCKET_NTERMS];
        bsim_rocket_xdot(
            P.data(), static_cast<int>(P.size()), get_doubles(prhs[2], BSIM_ROCKET_NX, "x"), u, nU,
            get_doubles(prhs[4], BSIM_ROCKET_AIR_N, "air"), get_scalar(prhs[5], "t"), xd, terms);
        plhs[0] = column(xd, BSIM_ROCKET_NX);
        plhs[1] = rocket_terms_struct(terms);
        plhs[2] = column(terms, BSIM_ROCKET_NTERMS);

    } else if (cmd == "isa_rocket_air") {
        check_args(nrhs, 4, cmd);
        double air[BSIM_ROCKET_AIR_N];
        bsim_isa_rocket_air(get_scalar(prhs[1], "h"), get_doubles(prhs[2], 3, "wind_i"),
                            get_doubles(prhs[3], 3, "gust_b"), air);
        plhs[0] = column(air, BSIM_ROCKET_AIR_N);

    } else if (cmd == "isa") {
        check_args(nrhs, 2, cmd);
        const mxArray *h = prhs[1];
        const double *hv = get_doubles(h, 0, "h");
        const std::size_t n = mxGetNumberOfElements(h);
        const mwSize nd = mxGetNumberOfDimensions(h);
        const mwSize *dims = mxGetDimensions(h);
        mxArray *out[4];
        double *o[4];
        for (int k = 0; k < 4; k++) {
            out[k] = mxCreateNumericArray(nd, dims, mxDOUBLE_CLASS, mxREAL);
            o[k] = mxGetDoubles(out[k]);
        }
        for (std::size_t i = 0; i < n; i++) {
            double r[4];
            bsim_isa(hv[i], r);
            for (int k = 0; k < 4; k++) {
                o[k][i] = r[k];
            }
        }
        for (int k = 0; k < 4; k++) {
            plhs[k] = out[k]; // matlab frees any that weren't asked for
        }

    } else if (cmd == "isa_air") {
        check_args(nrhs, 4, cmd);
        double air[8];
        bsim_isa_air(get_scalar(prhs[1], "h"), get_doubles(prhs[2], 3, "wind_i"),
                     get_doubles(prhs[3], 3, "gust_b"), air);
        plhs[0] = column(air, 8);

    } else if (cmd == "lookup1" || cmd == "lookup2") {
        const bool two = cmd == "lookup2";
        check_args(nrhs, two ? 6 : 4, cmd);
        const std::size_t nx = mxGetNumberOfElements(prhs[1]);
        const std::size_t ny = mxGetNumberOfElements(prhs[2]);
        const double *x = get_doubles(prhs[1], 0, "x");
        const double *y = get_doubles(prhs[2], two ? 0 : nx, "y");
        const double *z = two ? get_doubles(prhs[3], nx * ny, "Z (numel(x) by numel(y))") : nullptr;
        if (two) {
            check_table({x, nx}, {y, ny}, {z, nx * ny});
        } else {
            check_table({x, nx}, {y, ny});
        }

        const mxArray *xq = prhs[two ? 4 : 3];
        const double *xqv = get_doubles(xq, 0, "xq");
        const std::size_t n = mxGetNumberOfElements(xq);
        const double *yqv = two ? get_doubles(prhs[5], n, "yq (same size as xq)") : nullptr;

        mxArray *out = mxCreateNumericArray(mxGetNumberOfDimensions(xq), mxGetDimensions(xq),
                                            mxDOUBLE_CLASS, mxREAL);
        double *o = mxGetDoubles(out);
        for (std::size_t i = 0; i < n; i++) {
            o[i] = two ? bsim_lookup2(x, static_cast<int>(nx), y, static_cast<int>(ny), z, xqv[i],
                                      yqv[i])
                       : bsim_lookup1(x, y, static_cast<int>(nx), xqv[i]);
        }
        plhs[0] = out;

    } else {
        throw Error("unknown command '" + cmd + "', see the top of bsim_mex.cpp");
    }
}

} // namespace

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    static char message[1024];
    bool failed = false;
    try {
        dispatch(nlhs, plhs, nrhs, prhs);
    } catch (const std::exception &e) {
        std::snprintf(message, sizeof(message), "%s", e.what());
        failed = true;
    }
    if (failed) {
        mexErrMsgIdAndTxt("bsim:error", "%s", message);
    }
}
