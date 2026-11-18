#include "scenes/keyboard.hpp"
#include <algorithm>
#include <raylib.h>

namespace BSIM::scenes {

namespace {

// --- drone ---
Vector3d drone_keyboard_force(double m) {
    Vector3d f = Vector3d::Zero();

    double mult = 2.0 * m * 9.81;
    if (IsKeyDown(KEY_RIGHT_SHIFT)) {
        mult *= 1.2;
    }

    if (IsKeyDown(KEY_SPACE)) {
        f += Vector3d(0.0, 0.0, -mult);
    }
    if (IsKeyDown(KEY_LEFT_SHIFT)) {
        f += Vector3d(0.0, 0.0, mult);
    }

    return f;
}

Vector3d drone_keyboard_moment() {
    Vector3d T = Vector3d::Zero();

    double mult = 0.5;

    if (IsKeyDown(KEY_RIGHT_SHIFT)) {
        mult *= 1.2;
    }
    if (IsKeyDown(KEY_W)) {
        T.y() -= mult;
    }
    if (IsKeyDown(KEY_S)) {
        T.y() += mult;
    }
    if (IsKeyDown(KEY_A)) {
        T.x() -= mult;
    }
    if (IsKeyDown(KEY_D)) {
        T.x() += mult;
    }
    if (IsKeyDown(KEY_Q)) {
        T.z() -= mult;
    }
    if (IsKeyDown(KEY_E)) {
        T.z() += mult;
    }

    return T;
}

} // namespace

void read_keyboard(Drone &drone) {
    drone.input.F = drone_keyboard_force(drone.m);
    drone.input.M = drone_keyboard_moment();
}

// --- fixed wing ---
void read_keyboard(FixedWing &plane, double dt) {
    const double elevator = 0.1;      // rad on top of trim
    const double aileron = 0.1;       // rad on top of trim
    const double rudder = 0.15;       // rad on top of trim
    const double throttle_rate = 0.5; // per second

    FixedWingControls &controls = plane.controls;
    const FixedWingControls &trim = plane.trim;

    // positive elevator pitches the nose down (Cm_de < 0)
    controls.de = trim.de;
    if (IsKeyDown(KEY_W)) {
        controls.de += elevator;
    }
    if (IsKeyDown(KEY_S)) {
        controls.de -= elevator;
    }

    // positive aileron rolls right (Cl_da > 0)
    controls.da = trim.da;
    if (IsKeyDown(KEY_D)) {
        controls.da += aileron;
    }
    if (IsKeyDown(KEY_A)) {
        controls.da -= aileron;
    }

    // positive rudder yaws the nose left (Cn_dr < 0)
    controls.dr = trim.dr;
    if (IsKeyDown(KEY_Q)) {
        controls.dr += rudder;
    }
    if (IsKeyDown(KEY_E)) {
        controls.dr -= rudder;
    }

    // throttle stays where it's left
    if (IsKeyDown(KEY_UP)) {
        controls.dt += throttle_rate * dt;
    }
    if (IsKeyDown(KEY_DOWN)) {
        controls.dt -= throttle_rate * dt;
    }
    controls.dt = std::clamp(controls.dt, 0.0, 1.0);
}

} // namespace BSIM::scenes
