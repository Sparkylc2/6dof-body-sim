// src/main.cpp
#include "body.hpp"
#include "body_state.hpp"
#include "body_system.hpp"
#include "frames.hpp"
#include "render/render.hpp"
#include "world.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <raylib.h>

using namespace BSIM;

int main() {
    bool first_person = false;

    BSIM::rendering::init_rendering();
    BSIM::rendering::init_camera({-10.0, -6.0, -6.0}, {0, 0, 0}, {0, 0, -1}, 45.0);
    // the mesh is in raylib axes: x_rl is body y (span), y_rl is -body z, z_rl is -body x (chord)
    Mesh mesh = GenMeshCube(1.4224f, 0.05f, 0.3302f);
    // a small fin on top so the roll and which way is up are visible
    Mesh mesh_fin = GenMeshCube(0.02f, 0.15f, 0.2f);

    Model model = LoadModelFromMesh(mesh);
    Model model_fin = LoadModelFromMesh(mesh_fin);

    BSIM::World world;
    auto plane = std::make_unique<BSIM::FixedWing>(generate_sport_drone());
    //
    // auto plane = std::make_unique<BSIM::FixedWing>(generate_zagi_flying_wing());
    BSIM::FixedWing &body = *plane;
    world.bodies.bodies.push_back(std::move(plane));

    // level flight trim at 20 m/s, found numerically by solving fx = fz = m = 0 with theta = alpha
    // (for the zagi it was Va 18, alpha 0.071940, trim {-0.197314, 0, 0, 0.783290})
    const double Va_trim = 18, alpha_trim = 0.071940;
    // body.trim = {-0.197314, 0, 0, 0.783290};

    // starts 20 m up, flying north in trimmed level flight
    auto start = [&body, Va_trim, alpha_trim]() {
        body.state = BSIM::AircraftState();
        body.state.pos() = Vector3d(0.0, 0.0, -20.0);
        body.state.quat() = euler_to_quat(0.0, alpha_trim, 0.0);
        body.state.vel() =
            Vector3d(Va_trim * std::cos(alpha_trim), 0.0, Va_trim * std::sin(alpha_trim));
        body.controls = body.trim;
    };

    start();

    // physics steps at the fixed sensor rate Ts rather than frame time, so the analog
    // sensors get sampled at Ts. the frame time is fed in as whole steps with the left over
    // getting carried to the next frame
    //
    // todo: this should probably still be frame/fixed, but have a system for stepping the analog
    // sensors independently
    const double dt_sim = body.sensors.params.Ts;
    double t_behind = 0.0;

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_R)) {
            start();
        }
        if (IsKeyPressed(KEY_F)) {
            first_person = !first_person;
        }
        if (IsKeyPressed(KEY_T)) {
            world.atmosphere.turbulence = !world.atmosphere.turbulence;
        }
        const double dt = GetFrameTime();
        body.read_keyboard(dt);

        // ignore rendering until its caught up (if the sim falls behind the timestep)
        t_behind += std::min(dt, 0.25);
        while (t_behind >= dt_sim) {
            world.step(dt_sim);
            t_behind -= dt_sim;
        }
        const Eigen::Matrix3d R_b2i = rotation(body.state.quat());

        Matrix transform = rendering::to_rl_transform(R_b2i, body.state.pos());
        model.transform = transform;

        // fin sits slightly behind the middle of the wing and above it (-z in body)
        model_fin.transform =
            rendering::to_rl_transform(R_b2i, body.state.pos() + R_b2i * Vector3d(-0.1, 0, -0.075));

        if (first_person) {
            rendering::CAMERA.position = rendering::ned_to_rl(body.state.pos());
            rendering::CAMERA.target =
                rendering::ned_to_rl(body.state.pos() + R_b2i * Vector3d(1.0, 0.0, 0.0));
            rendering::CAMERA.up = rendering::ned_to_rl(R_b2i * Vector3d(0.0, 0.0, -1.0));
            rendering::CAMERA.fovy = 90.0f;
        } else {
            // chase camera behind and above, follows the heading but not the roll or pitch
            const double psi = quat_to_euler(body.state.quat()).z();
            const Vector3d behind(-std::cos(psi), -std::sin(psi), 0.0);
            const Vector3d pos = body.state.pos();
            BSIM::rendering::init_camera(pos + 6.0 * behind + Vector3d(0, 0, -1.5), pos, {0, 0, -1},
                                        45.0);
        }

        BSIM::rendering::begin();
        DrawModel(model, {0, 0, 0}, 1.0f, GRAY);
        DrawModel(model_fin, {0, 0, 0}, 1.0f, DARKGRAY);
        BSIM::rendering::end_3d();
        BSIM::rendering::draw_sensor_hud(body);
        BSIM::rendering::end();
    }
    CloseWindow();
}
