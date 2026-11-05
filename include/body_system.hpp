#pragma once
#include <Eigen/Dense>
#include <memory>
#include <vector>

#include "body.hpp"
#include "body_state.hpp"
#include "force_models.hpp"
#include "spacecraft.hpp"

namespace BSIM {
using Eigen::Ref;
using Eigen::Vector3d;
using Eigen::Vector4d;

// stores all the simulated bodies of one kind and does the physics updates
// B is the kind of body, it says which state it carries with B::State
template <class B> struct BodySystem {
    using State = typename B::State;
    using Derivative = StateVector<State::SIZE>;

    //  uses pointers so references to a body don't invalidate when more are added
    std::vector<std::unique_ptr<B>> bodies;

    // x_dot for body b at state x, each body type fills in its own stuff
    virtual Derivative derivatives(const B &b, const State &x) const = 0;

    // steps all bodies in the system forward, rk4 unless a system overrides it (only making it
    // virtual incase i want to switch integrators for different cases)
    virtual void step(double dt);

    // declaring the destructor would otherwise stop moves being generated
    // gotta love rule of 5
    BodySystem() = default;
    BodySystem(BodySystem &&) = default;
    BodySystem &operator=(BodySystem &&) = default;
    virtual ~BodySystem() = default;
};

// for symmetric aircraft flying through the air
struct UAVBodySystem : BodySystem<Body> {
    // gets the forces and moments from the body's wrench in the air it's in
    Derivative derivatives(const Body &b, const AircraftState &x) const override;

    // the equations of motion given the forces and moments (body frame)
    Derivative derivatives(const Body &b, const AircraftState &x, const Vector3d &F,
                           const Vector3d &T) const;
};

// for spacecraft (with no axes of symmetry, so full J tensor)
struct SpacecraftBodySystem : BodySystem<Spacecraft> {
    Derivative derivatives(const Spacecraft &s, const SpacecraftState &x) const override;
};

} // namespace BSIM
