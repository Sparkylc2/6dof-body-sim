// src/scenes/plane_scene.cpp
#include "autopilots/autopilot.hpp"
#include "autopilots/lateral_autopilot.hpp"
#include "body.hpp"
#include "body_state.hpp"
#include "body_system.hpp"
#include "frames.hpp"
#include "render/render.hpp"
#include "scenes/scenes.hpp"
#include "sensors.hpp" // wrap_angle
#include "trim.hpp"
#include "world.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <raylib.h>

namespace BSIM::scenes {

void run_plane_scene() {
    bool first_person = false;

    BSIM::rendering::init_camera({-10.0, -6.0, -6.0}, {0, 0, 0}, {0, 0, -1}, 45.0);
    Mesh mesh = GenMeshCube(1.4f, 0.05f, 0.3f);
    Mesh mesh_fin = GenMeshCube(0.02f, 0.15f, 0.2f);

    Model model = LoadModelFromMesh(mesh);
    Model model_fin = LoadModelFromMesh(mesh_fin);

    BSIM::World world;
    // auto plane = std::make_unique<BSIM::FixedWing>(generate_zagi_flying_wing());
    //
    auto plane = std::make_unique<BSIM::FixedWing>(generate_sport_drone());
    BSIM::FixedWing &body = *plane;
    world.bodies.bodies.push_back(std::move(plane));

    // level flight trim at 18 m/s
    const double Va_trim = 18;
    const BSIM::LongitudinalTrim trim = BSIM::trim_longitudinal(body, Va_trim);
    if (!trim.converged) {
        TraceLog(LOG_WARNING, "trim didn't converge, residual %g", trim.residual);
    }
    body.trim = trim.controls;

    auto start = [&body, &trim]() {
        body.state = trim.state(Vector3d(0.0, 0.0, -20.0));
        body.controls = body.trim;
    };

    start();

    // lateral autopilot
    BSIM::LateralDesignParams lat;
    lat.da_max = 45.0 * DEG;    // seems an appropriate angle
    lat.e_phi_max = 30.0 * DEG; // same here
    lat.zeta_phi = 0.9;         // this shoulld be well damped
    lat.ki_phi = 0.1;   // a small value should just trim steady-state and not destabilise anything
    lat.W_chi = 5.0;    // textbook says this should be good
    lat.zeta_chi = 1.0; // larger is better since this gives more bandwidth/less overshoot
    lat.Vg = Va_trim;   // trim condition
    lat.phi_max = 30.0 * DEG;    // reasonable bank limit
    lat.dr_max = 30.0 * DEG;     // reasonable rudder limit
    lat.e_beta_max = 10.0 * DEG; // i'm guessing here i'll change if its shit
    lat.zeta_beta = 0.9;         // should also be well damped

    BSIM::LateralAutopilot lat_ap;
    lat_ap.clock.Ts = body.sensors.params.Ts;
    lat_ap.design(body.params, body.m, body.J, BSIM::RHO_DEFAULT, Va_trim, lat);

    bool ap_on = false;
    double chi_c = 0.0;
    BSIM::FixedWingControls ap_u; // last autopilot outputs
    double t_sim = 0.0;

    // physics steps at the fixed sensor rate Ts rather than frame time, so the analog
    // sensors get sampled at Ts. the frame time is fed in as whole steps with the left over
    // getting carried to the next frame
    //
    // todo: decouple frame time from the sensor Ts
    const double dt_sim = body.sensors.params.Ts;
    double t_behind = 0.0;

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_R)) {
            start();
            ap_on = false;
        }
        if (IsKeyPressed(KEY_P)) {
            ap_on = !ap_on;
            if (ap_on) {
                chi_c = BSIM::true_control_state(body.state, body.air).chi;
                lat_ap.reset();
                lat_ap.clock.next = t_sim;
                ap_u = body.controls;
            }
        }
        if (ap_on && IsKeyPressed(KEY_LEFT)) {
            chi_c = BSIM::wrap_angle(chi_c - 15.0 * DEG);
        }
        if (ap_on && IsKeyPressed(KEY_RIGHT)) {
            chi_c = BSIM::wrap_angle(chi_c + 15.0 * DEG);
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
            if (ap_on) {
                // flies off the true state for now
                const BSIM::ControlState s = BSIM::true_control_state(body.state, body.air);
                lat_ap.step(t_sim, chi_c, s.phi, s.chi, s.beta, s.p, ap_u);
                body.controls.da = ap_u.da;
                if (lat_ap.has_rudder) {
                    body.controls.dr = ap_u.dr;
                }
            }
            world.step(dt_sim);
            t_sim += dt_sim;
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
        if (ap_on) {
            DrawText(TextFormat("lateral ap ON  chi_c %6.1f deg     phi_c %6.1f deg     [p]"
                                "[<- or ->] for course",
                                chi_c / DEG, lat_ap.phi_c / DEG),
                     10, GetScreenHeight() - 30, 16, DARKGREEN);
        } else {
            DrawText("lateral ap off  [p]", 10, GetScreenHeight() - 30, 16, DARKGRAY);
        }
        BSIM::rendering::end();
    }
    UnloadModel(model);
    UnloadModel(model_fin);
}

} // namespace BSIM::scenes
