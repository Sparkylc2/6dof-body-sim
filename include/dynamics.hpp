#pragma once
#include "body_state.hpp"
#include "force_models/rocket_model.hpp"
#include "force_models/spacecraft_model.hpp"
#include <Eigen/Dense>

// the equations of motion as plain functions, everything they need comes in as arguments and
// nothing is kept between calls (the BodySystems call these, and so can anything outside the sim
// like the matlab bridge)
namespace BSIM {
using Eigen::Vector3d;
using Eigen::Vector4d;

// --- aircraft ---
// x_dot for a symmetric aircraft given the forces F and moments T (both body frame)
// m is the mass and J is [Jx, Jy, Jz, Jxz]
// from eqn B.2 and B.4
StateVector<AircraftState::SIZE> aircraft_derivatives(double m, const Vector4d &J,
                                                      const AircraftState &x, const Vector3d &F,
                                                      const Vector3d &T);

// --- spacecraft ---
// x_dot for a spacecraft with reaction wheels, thrusters and a fuel tank
// (works out its own forces and moments from the params, environment and controls)
StateVector<SpacecraftState::SIZE> spacecraft_derivatives(const SpacecraftParams &p,
                                                          const SpacecraftEnvironment &env,
                                                          const SpacecraftControls &u,
                                                          const SpacecraftState &x);

// --- rocket ---
// x_dot for a variable mass rocket with a gimballed engine and grid fins
// (works out its own forces and moments, see force_models/rocket_model.hpp)
//   u is where the actuators actually are, air comes from the atmosphere / wind models, and t is
//   the sim time (only thrust curve engines use it, against u.t_ignition)
//   terms, if given, gets every force and moment it used (for logging and the accelerometer)
StateVector<RocketState::SIZE> rocket_derivatives(const RocketParams &p,
                                                  const RocketEnvironment &env,
                                                  const RocketControls &u, const RocketAir &air,
                                                  const RocketState &x, double t,
                                                  RocketTerms *terms = nullptr);

} // namespace BSIM
