/*
 * bsim_c.h: the c interface to the sim's models, for simulink (C Caller / C Function blocks)
 * and anything else that can call c
 *
 * every function is stateless (same inputs, same outputs, nothing kept between calls), so time,
 * integration, gusts, actuators, sensors and control all stay on the calling side
 *
 * everything is plain double arrays, matrices are column major (same as matlab)
 * the parameter vectors P come packed, get them from matlab with bsim.fixedwing.params /
 * bsim.spacecraft.params (or the factory functions below)
 *
 * layouts:
 *   aircraft state x [13]:   [pn pe pd | u v w | e0 e1 e2 e3 | p q r]   (u v w = ground vel, body)
 *   fixed wing controls [4]: [de da dr dt]   (where the surfaces actually are, rad, throttle 0-1)
 *   air [8]:                 [rho a | wind_i (3, inertial) | gust_b (3, body)]
 *   spacecraft state x [17]: [r (3) | v (3) | q (4) | w (3) | h_wheels (3) | m]   (r, v inertial)
 *   rocket state x [14]:     [pn pe pd | vn ve vd | e0 e1 e2 e3 | p q r | m]  (cm pos/vel, NED)
 *   rocket controls [9 + n_fins], or [5 + n_fins] without the gimbal's motion (taken as still):
 *                            [gimbal_p gimbal_y throttle t_ignition fins_deployed | fin_1..fin_n |
 *                             gimbal_p_dot gimbal_y_dot gimbal_p_ddot gimbal_y_ddot]
 *                            (where the actuators actually are, rad, throttle 0-1, t_ignition is
 *                             the sim time a thrust curve motor was lit, inf for not yet. the
 *                             gimbal's rates and accelerations, from its servos, drive the swinging
 *                             engine's dynamics, see force_models/tvc_model.hpp)
 *   rocket air [9]:          [rho a p | wind_i (3, NED) | gust_b (3, body)]
 */
#ifndef BSIM_C_H
#define BSIM_C_H

#ifdef __cplusplus
extern "C" {
#endif

#define BSIM_AIRCRAFT_NX 13
#define BSIM_FIXEDWING_NU 4
#define BSIM_AIR_N 8
#define BSIM_SPACECRAFT_NX 17

/* [m, Jx, Jy, Jz, Jxz, then every FixedWingParams field in order] */
#define BSIM_FIXEDWING_NP 47

/* [m_empty, J_empty (9), r_tank (3), R_tank, wheel A (9), J_w, tau_max, h_max,
 *  env gravity, env j2, env gravity_gradient (1 = on), n_thrusters,
 *  then per thruster: pos (3), dir (3), T_max, Isp]                                  */
#define BSIM_SPACECRAFT_NP_BASE 30
#define BSIM_THRUSTER_NP 8
#define BSIM_SPACECRAFT_NP(n_thrusters) (BSIM_SPACECRAFT_NP_BASE + BSIM_THRUSTER_NP * (n_thrusters))

#define BSIM_ROCKET_NX 14
#define BSIM_ROCKET_NU(n_fins) (9 + (n_fins))
#define BSIM_ROCKET_NU_STILL(n_fins) (5 + (n_fins)) /* the gimbal's motion left out */
#define BSIM_ROCKET_AIR_N 9

/* the rocket's P is variable length, a fixed head then the lists it says it has:
 *   [0]  m_dry x_dry Ix_dry It_dry
 *   [4]  d A_base A_motor C_dc CD_flat
 *   [9]  engine: kind (0 throttleable, 1 thrust curve) c_eff p_test T_vac_max Isp_vac
 *                throttle_min A_e x_G L_n m_n l_n I_na I_nt
 *   [22] fins: n x r phi0 S CN_alpha CD0
 *   [29] env: gravity aero jet_damping J_dot coriolis (1 = on)
 *   [34] counts: n_propellant n_components n_slices n_thrust_curve n_friction n_pressure
 *                n_parasitic
 *   [41] then in that order:
 *        propellant  6 each [kind (0 liquid, 1 solid) share x R rho L]
 *        components  2 each [CN_alpha x_cp]
 *        slices      3 each [x d dx]
 *        each table (thrust curve, friction, pressure, parasitic): all its x, then all its y
 * an empty table is fine (it counts as zero)                                              */
#define BSIM_ROCKET_NP_HEAD 41

/* what came out of the rocket's x_dot (body frame unless it says otherwise) */
#define BSIM_ROCKET_NTERMS 41
/*   [0]  m x_cm Ix It           mass properties right now (x_cm is a station, see rocket_model.hpp)
 *   [4]  T m_dot_e              thrust [N], exhaust mass flow [kg/s]
 *   [6]  V mach qbar alpha_t    air data (alpha_t is the total angle of attack, 0 to pi)
 *   [10] F_body (3)             every non gravitational force, F_body / m is what an
 *                               accelerometer at the cm reads
 *   [13] M_body (3)             every moment about the cm
 *   [16] F_aero (3), [19] M_aero (3)
 *   [22] g_i (3)                gravity (NED)
 *   [25] n_T (3)                direction the thrust pushes (follows the gimbal, even when off)
 *   [28] rho_e (3)              nozzle exit from the cm
 *   [31] v_r (3)                velocity relative to the air
 *   [34] M_gimbal (3)           every moment from the swinging engine (tail wags dog, coriolis,
 *                               changing geometry, swinging exit), already part of M_body
 *   [37] hinge (2)              torque the gimbal servos have to supply [pitch, yaw] [N m]
 *   [39] nozzle_load (2)        the vehicle's load on the nozzle (what the TVC's load input
 *                               wants, a step late is fine) [N m]                               */

/* ---- aircraft ---- */

/* x_dot for a symmetric aircraft given body frame forces F and moments M, J = [Jx Jy Jz Jxz] */
void bsim_aircraft_xdot(double m, const double J[4], const double x[13], const double F[3],
                        const double M[3], double xdot[13]);

/* ---- fixed wing ---- */

/* forces and moments (body frame) from gravity, aero and propulsion */
void bsim_fixedwing_wrench(const double P[47], const double x[13], const double u[4],
                           const double air[8], double F[3], double M[3]);

/* x_dot, plus the forces and moments it came from (the accelerometer model wants them) */
void bsim_fixedwing_xdot(const double P[47], const double x[13], const double u[4],
                         const double air[8], double xdot[13], double F[3], double M[3]);

/* d(x_dot)/dx (13 x 13) and d(x_dot)/du (13 x 4), central finite differences */
void bsim_fixedwing_jac(const double P[47], const double x[13], const double u[4],
                        const double air[8], double A[169], double B[52]);

/* the packed parameters of the built in planes */
void bsim_fixedwing_zagi(double P[47]);
void bsim_fixedwing_sport(double P[47]);

/* ---- spacecraft ---- */

/* x_dot, throttle has one entry per thruster in P
 * returns 0, or -1 if nP doesn't match the thruster count inside P (xdot is left alone then) */
int bsim_spacecraft_xdot(const double *P, int nP, const double x[17], const double *throttle,
                         const double wheel_torque[3], double xdot[17]);

/* the packed parameters and starting state of the built in cubesat
 * returns how many entries P needs (nothing is written if nP_max is too small) */
int bsim_spacecraft_cubesat(double *P, int nP_max, double x0[17]);

/* ---- rocket ---- */

/* x_dot, and the terms it came from (see BSIM_ROCKET_NTERMS)
 * returns 0, or -1 if P isn't a valid packed rocket or nU isn't 9 (or 5) + its fin count
 * (xdot and terms are filled with NaN then, so a sim using them goes obviously wrong) */
int bsim_rocket_xdot(const double *P, int nP, const double x[14], const double *u, int nU,
                     const double air[9], double t, double xdot[14], double terms[41]);

/* the packed parameters and starting state (on the pad, nose up, full) of the example rocket
 * returns how many entries P needs (nothing is written if nP_max is too small) */
int bsim_rocket_example(double *P, int nP_max, double x0[14]);

/* the rocket's air vector for altitude h from ISA, with the given winds */
void bsim_isa_rocket_air(double h, const double wind_i[3], const double gust_b[3], double air[9]);

/* ---- hydraulic TVC (force_models/tvc_model.hpp, the Ember TVC actuator notes) ----
 * four hydraulic pistons swing the nozzle. it keeps no state: the 12 states live with the caller
 * (a Unit Delay in simulink), and bsim_tvc_step moves them on by dt, sub stepping the stiff
 * hydraulics at its own small step inside                                                     */
#define BSIM_TVC_NX 12   /* [delta (2) | delta_dot (2) | P_L (4) | u_v (4)] */
#define BSIM_TVC_NP 23   /* every TvcParams field in order, see BSIM_TVC_FIELDS in bridge.hpp */
#define BSIM_TVC_NOUT 18
/*   [0]  delta (2)       where the nozzle is [rad]
 *   [2]  delta_dot (2)   [rad/s]
 *   [4]  delta_ddot (2)  [rad/s^2] (from bsim_tvc_step: the average over the step, which is what
 *                        the rocket's tail wags dog term wants)
 *   [6]  l (4)           piston lengths [m]
 *   [10] l_cmd (4)       what each piston is being asked for [m]
 *   [14] F (4)           piston forces [N]                                                   */

/* the placeholder Ember TVC, with the example rocket's engine as the nozzle */
void bsim_tvc_ember(double P[23]);

/* moves the TVC on by dt with the inputs held: the gimbal command [rad], the vehicle's load on the
 * nozzle [N m] (the rocket terms' nozzle_load) and the thrust [N] (for the thrust offset) */
void bsim_tvc_step(const double P[23], const double s0[12], const double delta_cmd[2],
                   const double load[2], double thrust, double dt, double s1[12], double out[18]);

/* the 12 state derivatives (for analysis, simulations use bsim_tvc_step) */
void bsim_tvc_xdot(const double P[23], const double s[12], const double delta_cmd[2],
                   const double load[2], double thrust, double sdot[12], double out[18]);

/* piston lengths [m] and the jacobian l_dot = J delta_dot (4 x 2, column major) at delta */
void bsim_tvc_kinematics(const double P[23], const double delta[2], double l[4], double J[8]);

/* gimbal angles from four measured piston lengths (gauss newton from delta0),
 * returns the leftover length misfit [m] (big means a bad sensor) */
double bsim_tvc_forward_kinematics(const double P[23], const double l[4], const double delta0[2],
                                   double delta[2]);

/* ---- atmosphere ---- */

/* ISA at geometric altitude h [m], out = [T p rho a] */
void bsim_isa(double h, double out[4]);

/* the air vector for altitude h from ISA, with the given winds */
void bsim_isa_air(double h, const double wind_i[3], const double gust_b[3], double air[8]);

/* ---- lookup tables ---- */
/* linear, held flat past the ends, breakpoints strictly increasing (not checked here)
 * z is nx by ny column major (z[i + nx*j] is at (x[i], y[j]))                          */
double bsim_lookup1(const double *x, const double *y, int n, double xq);
double bsim_lookup2(const double *x, int nx, const double *y, int ny, const double *z, double xq,
                    double yq);

#ifdef __cplusplus
}
#endif

#endif /* BSIM_C_H */
