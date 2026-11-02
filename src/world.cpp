#include "world.hpp"
#include "frames.hpp"

namespace uav {

void World::step(double dt) {
    for (auto &bp : bodies.bodies) {
        Body &b = *bp;

        // airspeed from the current state for the gust filters
        // you need to rotate the steady wind into body axes cuz its in inertial rn
        const Vector3d wind_b = rotation(b.state.quat()).transpose() * atmosphere.wind_i;
        // Vg = Va + Vw => Va = Vg - Vw, as state velocity is ground speed
        const double Va = (b.state.vel() - wind_b).norm();

        // update the bodies air state to match the current atmosphere
        b.air.rho = atmosphere.rho;       // for aircraft we'd prob add a function here using alt
        b.air.wind_i = atmosphere.wind_i; // same here
        // gusts advance once per step and are then held through all the rk4 stages
        b.air.gust_b = atmosphere.turbulence && b.gust ? b.gust->step(Va, dt, atmosphere.rng)
                                                       : Vector3d::Zero();
    }

    bodies.step(dt);
    t += dt;

    // the sensors get the state after integration/end of the physics step
    for (auto &bp : bodies.bodies) {
        bp->sensors.update(*bp, t, sensor_rng);
    }
}
} // namespace uav
