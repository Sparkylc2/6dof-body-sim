#pragma once
#include "body_state.hpp"
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

} // namespace BSIM
