#pragma once
#include "air_system.hpp"
#include "body_system.hpp"

namespace uav {

// describes the world the bodies live in
struct World {
    AirSystem atmosphere; // atmospheric model for wind and turbulence
    UAVBodySystem bodies; // dynamics model for the bodies
    double t = 0.0;

    // separate from the atmosphere's so changing the sensors doesn't change the turbulence
    std::mt19937 sensor_rng{1};

    // steps the atmosphere and the bodies, then samples each body's sensors
    void step(double dt);
};
} // namespace uav
