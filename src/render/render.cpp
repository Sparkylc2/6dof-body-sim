#include "render/render.hpp"
#include "frames.hpp"
#include "sensors.hpp"

namespace BSIM::rendering {

Camera3D CAMERA = {};

void init_rendering() {
    InitWindow(1280, 720, "bsim");
    SetTargetFPS(60);
}

void init_camera(const Eigen::Vector3d &pos, const Eigen::Vector3d &target,
                 const Eigen::Vector3d &up, double fovy) {
    CAMERA.position = ned_to_rl(pos);
    CAMERA.target = ned_to_rl(target);
    CAMERA.up = ned_to_rl(up);
    CAMERA.fovy = static_cast<float>(fovy);
    CAMERA.projection = CAMERA_PERSPECTIVE;
}

void begin(Color background, bool ground) {
    BeginDrawing();
    ClearBackground(background);
    BeginMode3D(CAMERA);
    if (ground) {
        DrawGrid(200, 5.0f); // 1 km across so there's still ground to see at flying speed
    }
}

void end_3d() { EndMode3D(); }

void end() {
    DrawFPS(10, 10);
    EndDrawing();
}

void draw_sensor_hud(const Body &b) {
    const AircraftState &x = b.state;
    const SensorReadings &y = b.sensors.y;

    // the noise free versions of each reading
    const Vector3d f = specific_force(b.wrench(x, b.air), b.m, x.quat());
    const AirData air = air_data(x.vel(), x.quat(), b.air);
    const Vector2d Vg_chi = ground_speed_course(rotation(x.quat()) * x.vel());
    const double psi = quat_to_euler(x.quat()).z();
    const double r2d = 1.0 / DEG;

    int row = 40;
    auto line = [&row](const char *text) {
        DrawText(text, 10, row, 16, DARKGRAY);
        row += 20;
    };

    line("sensor          measured                  (true)");
    line(TextFormat("gyro deg/s     %7.2f %7.2f %7.2f   (%7.2f %7.2f %7.2f)", y.gyro.x() * r2d,
                    y.gyro.y() * r2d, y.gyro.z() * r2d, x.rates().x() * r2d, x.rates().y() * r2d,
                    x.rates().z() * r2d));
    line(TextFormat("accel m/s^2    %7.2f %7.2f %7.2f   (%7.2f %7.2f %7.2f)", y.accel.x(),
                    y.accel.y(), y.accel.z(), f.x(), f.y(), f.z()));
    line(TextFormat("abs pres Pa    %7.1f   (%7.1f)", y.abs_pres,
                    abs_pressure(b.air.rho, -x.pos().z())));
    line(TextFormat("diff pres Pa   %7.1f   (%7.1f)", y.diff_pres, diff_pressure(air.rho, air.Va)));
    line(TextFormat("compass deg    %7.2f   (%7.2f)", y.compass * r2d, psi * r2d));
    line(TextFormat("gps n e h m    %7.1f %7.1f %7.1f   (%7.1f %7.1f %7.1f)", y.gps_pos.x(),
                    y.gps_pos.y(), y.gps_pos.z(), x.pos().x(), x.pos().y(), -x.pos().z()));
    line(TextFormat("gps Vg m/s     %7.2f   (%7.2f)", y.gps_Vg, Vg_chi.x()));
    line(TextFormat("gps course deg %7.2f   (%7.2f)", y.gps_course * r2d, Vg_chi.y() * r2d));
}

} // namespace BSIM::rendering
