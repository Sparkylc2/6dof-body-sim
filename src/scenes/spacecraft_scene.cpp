// src/scenes/spacecraft_scene.cpp
#include "autopilots/attitude_control.hpp"
#include "body_system.hpp"
#include "frames.hpp"
#include "render/render.hpp"
#include "scenes/scenes.hpp"
#include "sensors.hpp" // DEG
#include "spacecraft.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <random>
#include <raylib.h>

namespace BSIM::scenes {

void run_spacecraft_scene() {
    SpacecraftBodySystem sys;
    sys.bodies.push_back(std::make_unique<Spacecraft>(generate_cubesat()));
    Spacecraft &sc = *sys.bodies[0];
    const Spacecraft initial = sc;

    AttitudeController ctrl;
    ctrl.clock.Ts = 0.1;
    ctrl.design(sc.params, sc.state.mass(), AttitudeDesignParams{});

    // points the body z axis (b) at a random inertial direction
    std::mt19937 rng(std::random_device{}());
    std::normal_distribution<double> n(0.0, 1.0);
    auto random_direction = [&rng, &n]() { return Vector3d(n(rng), n(rng), n(rng)).normalized(); };

    AttitudeCommand cmd;
    cmd.mode = AttitudeCommand::Mode::Point;
    cmd.b = Vector3d::UnitZ();
    cmd.d_i = random_direction();

    const double dt_sim = 0.05;
    double t_sim = 0.0;
    double t_behind = 0.0;

    // drawn 5x so it's visible against the
    const double S = 5.0;
    const Vector3d centre(0.0, 0.0, -3.0);
    const double line_length = 0.4 * S;

    // cube inputs are in raylib's axes, which are body y, -z, -x once transformed
    Model box = LoadModelFromMesh(GenMeshCube(S * 0.226f, S * 0.366f, S * 0.10f));
    rendering::init_camera(centre + Vector3d(-7.0, -4.0, -2.5), centre, {0, 0, -1}, 45.0);

    // a line from the centre along a unit direction (in inertial axes)
    auto draw_direction = [&](const Vector3d &dir, Color c) {
        const Vector3 a = rendering::ned_to_rl(centre);
        const Vector3 b = rendering::ned_to_rl(centre + line_length * dir);
        DrawCylinderEx(a, b, 0.03f, 0.03f, 8, c);
        DrawSphere(b, 0.1f, c);
    };

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_SPACE)) {
            cmd.d_i = random_direction();
        }
        if (IsKeyPressed(KEY_R)) {
            sc = initial;
            t_sim = 0.0;
            ctrl.reset();
            ctrl.clock.next = 0.0;
            cmd.d_i = random_direction();
        }

        // same as the plane, fixed physics steps with the leftover carried to the next frame
        t_behind += std::min(static_cast<double>(GetFrameTime()), 0.25);
        while (t_behind >= dt_sim) {
            ctrl.step(t_sim, cmd, sc.state.quat(), sc.state.rates(), sc.controls);
            sys.step(dt_sim);
            t_sim += dt_sim;
            t_behind -= dt_sim;
        }

        const Matrix3d R = rotation(sc.state.quat());
        const Vector3d b_i = R * cmd.b;
        const double err = std::acos(std::clamp(b_i.dot(cmd.d_i), -1.0, 1.0));

        rendering::begin();
        box.transform = rendering::to_rl_transform(R, centre);
        DrawModel(box, {0, 0, 0}, 1.0f, LIGHTGRAY);
        DrawModelWires(box, {0, 0, 0}, 1.0f, DARKGRAY);
        draw_direction(b_i, BLUE);       // b, where it's pointing
        draw_direction(cmd.d_i, ORANGE); // where it should be pointing
        rendering::end_3d();

        DrawText(TextFormat("error %6.2f deg     blue b, orange target     [space] new target  "
                            "[r] reset",
                            err / DEG),
                 10, GetScreenHeight() - 30, 16, DARKGRAY);
        rendering::end();
    }
    UnloadModel(box);
}

} // namespace BSIM::scenes
