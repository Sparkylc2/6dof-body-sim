#pragma once
#include <limits>

namespace BSIM {

// --- actuator model ---
// servo that takes a commanded deflection and moves the surface towards it
//   ideal:        pos = cmd
//   first order:  pos' = (cmd - pos) / tau
//   second order: pos'' = wn^2 (cmd - pos) - 2 zeta wn pos'
// on top of the linear dynamics the rate and position are clamped (making it nonlinear)
struct Actuator {
    enum class Kind { Ideal, FirstOrder, SecondOrder };
    Kind kind = Kind::Ideal;

    double tau = 0.05;            // first order time constant [s]
    double wn = 30.0, zeta = 0.7; // second order natural freq [rad/s] and damping
    double rate_max = std::numeric_limits<double>::infinity(); // [rad/s]
    double pos_max = std::numeric_limits<double>::infinity();  // [rad]

    double pos = 0.0;  // where the surface actually is
    double rate = 0.0; // how fast it's moving
    double acc = 0.0;  // its average acceleration over the last step (the rocket's gimbal dynamics
                       // use it for the tail wags dog reaction)

    // holds cmd over dt and moves the surface, returns the new pos
    // (sub steps internally, gotta love stiffness)
    double step(double cmd, double dt);

    // puts the surface at pos, at rest
    void reset(double p = 0.0);

  private:
    double step_motion(double cmd, double dt); // the step itself, step() adds acc around it
  public:

    bool ideal() const { return kind == Kind::Ideal; }
};

} // namespace BSIM
