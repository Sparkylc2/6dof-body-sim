#pragma once
// the c++ side of the bridge, shared by the c api (bsim_c.cpp) and the mex (bsim_mex.cpp)
// this is the only place the packed layouts in bsim_c.h are written down in code
#include "actuators.hpp"
#include "air_system.hpp"
#include "body_state.hpp"
#include "force_models/fixed_wing_model.hpp"
#include "force_models/rocket_model.hpp"
#include "force_models/spacecraft_model.hpp"
#include "force_models/tvc_model.hpp"
#include <Eigen/Dense>
#include <vector>

namespace BSIM::Bridge {

// every FixedWingParams field, in the order they're packed (after [m, Jx, Jy, Jz, Jxz])
// (adding a field to FixedWingParams without adding it here fails a static_assert in bsim_c.cpp)
// clang-format off
#define BSIM_FIXEDWING_PARAM_FIELDS(X)                                                             \
    X(S) X(b) X(c) X(e)                                                                            \
    X(CL0) X(CL_alpha) X(CL_q) X(CL_de)                                                            \
    X(CD0) X(CD_alpha) X(CD_q) X(CD_de) X(CD_p)                                                    \
    X(Cm0) X(Cm_alpha) X(Cm_q) X(Cm_de)                                                            \
    X(M) X(alpha0)                                                                                 \
    X(CY0) X(CY_beta) X(CY_p) X(CY_r) X(CY_da) X(CY_dr)                                            \
    X(Cl0) X(Cl_beta) X(Cl_p) X(Cl_r) X(Cl_da) X(Cl_dr)                                            \
    X(Cn0) X(Cn_beta) X(Cn_p) X(Cn_r) X(Cn_da) X(Cn_dr)                                            \
    X(S_prop) X(C_prop) X(k_motor) X(k_Tp) X(k_Omega)
// clang-format on

constexpr int FIXEDWING_HEAD = 5; // m, Jx, Jy, Jz, Jxz

// everything the fixed wing functions need (a FixedWing without any of the state)
struct FixedWingDef {
    double m = 0.0;
    Eigen::Vector4d J = Eigen::Vector4d::Zero();
    FixedWingParams params;
};

FixedWingDef unpack_fixedwing(const double *P);
void pack_fixedwing(const FixedWingDef &d, double *P);

FixedWingDef zagi();
FixedWingDef sport();

struct SpacecraftDef {
    SpacecraftParams params;
    SpacecraftEnvironment env;
};

int spacecraft_np(int n_thrusters);
// false if nP doesn't match the thruster count stored in P
bool unpack_spacecraft(const double *P, int nP, SpacecraftDef &out);
// P needs spacecraft_np(thrusters) entries
void pack_spacecraft(const SpacecraftDef &d, double *P);

// everything the rocket functions need apart from the controls
struct RocketDef {
    RocketParams params;
    RocketEnvironment env;
};

// how long the packed P is for this rocket
int rocket_np(const RocketDef &d);
// throws std::invalid_argument saying what's wrong if P isn't a valid packed rocket
RocketDef unpack_rocket(const double *P, int nP);
std::vector<double> pack_rocket(const RocketDef &d);

// the actuators that go with the example rocket
struct RocketExample {
    RocketDef def;
    RocketState x0;
    Actuator fin;  // the grid fin servos
    TvcParams tvc; // the hydraulic TVC swinging the engine
};
RocketExample rocket_example();

// with_motion: u carries the gimbal's rates and accelerations after the fins (see bsim_c.h)
RocketControls unpack_rocket_controls(const double *u, int n_fins, bool with_motion);
RocketAir unpack_rocket_air(const double *air);

AirState unpack_air(const double *air);

// every TvcParams field, in the order they're packed
// (adding a field to TvcParams without adding it here fails a static_assert in bsim_c.cpp)
// clang-format off
#define BSIM_TVC_FIELDS(X)                                                                         \
    X(x_a) X(r_a) X(x_b) X(r_b) X(phi0) X(stroke)                                                  \
    X(A_p) X(P_s) X(Q_max) X(beta) X(V_t) X(C_tp) X(tau_v) X(K_p)                                  \
    X(m_n) X(l_n) X(I_na) X(I_nt)                                                                  \
    X(offset_p) X(offset_y) X(friction) X(damping)                                                 \
    X(h)
// clang-format on

TvcParams unpack_tvc(const double *P);
void pack_tvc(const TvcParams &p, double *P);

} // namespace BSIM::Bridge
