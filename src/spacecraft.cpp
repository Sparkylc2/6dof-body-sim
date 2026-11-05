#include "spacecraft.hpp"
#include <numbers>

namespace BSIM {

Spacecraft generate_cubesat() {
    Spacecraft s;
    SpacecraftParams &p = s.params;

    // empty box with dimensions
    // 0.10 x 0.226 x 0.366 m,
    // m/12 (b^2 + c^2) on the diagonal,
    // plus a bit of inertia on the diagonals (just to see something not symmetric)
    p.m_empty = 11.0;
    p.J_empty << 0.170, -0.001, 0.002, //
        -0.001, 0.130, -0.0005,        //
        0.002, -0.0005, 0.056;

    // tank is on the z axis below the empty cm, so the cm will stay on the thrust line as it gets
    // emptied
    p.r_tank = Vector3d(0.0, 0.0, -0.08);
    p.R_tank = 0.04;

    // small cubesat wheels (quick search says ~10000 rpm at l_max is reasonable)
    p.wheels.A = Matrix3d::Identity();
    p.wheels.J_w = 2.0e-5;
    p.wheels.tau_max = 2.0e-3;
    p.wheels.h_max = 0.02;

    // one main thruster on the -z face pushing along +z
    p.thrusters.push_back(
        {.pos = Vector3d(0.0, 0.0, -0.183), .dir = Vector3d::UnitZ(), .T_max = 1.0, .Isp = 220.0});
    s.controls.throttle.assign(p.thrusters.size(), 0.0);

    const auto [r, v] = circular_orbit(500.0e3, 51.6 * std::numbers::pi / 180.0);
    s.state.pos() = r;
    s.state.vel() = v;
    s.state.mass() = p.m_empty + 1.0; // ssuming a full tank

    return s;
}

} // namespace BSIM
