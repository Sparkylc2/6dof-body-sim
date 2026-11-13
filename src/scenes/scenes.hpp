#pragma once

namespace BSIM::scenes {

// each scene runs its own loop until the window is closed
// (the window is opened and closed by main)

// plane in level flight over the grid
void run_plane_scene();

// the cubesat (drawn as a box) using the reaction wheel attitude controller to point its body z
// axis at a random direction
void run_spacecraft_scene();

// initially just runs the different trajectories to develop the
// SINDy model on phugoid data
void run_sindy_scene();

} // namespace BSIM::scenes
