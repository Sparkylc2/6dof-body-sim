#pragma once
#include "craft.hpp"

// keyboard controls for the bodies, kept out of bsim_core so the physics doesn't need raylib
// (call once per frame before stepping)
namespace BSIM::scenes {

// space/left shift thrust up/down, w/s pitch, a/d roll, q/e yaw, right shift boosts everything
// (fills drone.input)
void read_keyboard(Drone &drone);

// w/s elevator (nose down/up), a/d aileron (roll left/right), q/e rudder (yaw left/right),
// up/down arrows throttle
void read_keyboard(FixedWing &plane, double dt);

} // namespace BSIM::scenes
