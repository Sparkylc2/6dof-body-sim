#pragma once
#include "body.hpp"
#include <Eigen/Dense>
#include <raylib.h>
#include <raymath.h>

namespace BSIM::rendering {

extern Camera3D CAMERA; // global camera

// raylib uses a different coord convention (y is up)
inline Vector3 ned_to_rl(const Eigen::Vector3d &p) {
    return {static_cast<float>(p.y()), static_cast<float>(-p.z()), static_cast<float>(-p.x())};
}

// converts from the Eigen NED frame to raylibs y up frame, applying the rotation and translation
inline Matrix to_rl_transform(const Eigen::Matrix3d &R, const Eigen::Vector3d &p) {
    Eigen::Matrix3d P;
    P << 0, 1, 0, // x_rl = e
        0, 0, -1, // y_rl = −d
        -1, 0, 0; // z_rl = −n
    Eigen::Matrix3d R_rl = P * R * P.transpose();
    Eigen::Vector3d p_rl = P * p;

    ::Matrix m = MatrixIdentity();
    m.m0 = R_rl(0, 0);
    m.m4 = R_rl(0, 1);
    m.m8 = R_rl(0, 2);
    m.m12 = p_rl.x();
    m.m1 = R_rl(1, 0);
    m.m5 = R_rl(1, 1);
    m.m9 = R_rl(1, 2);
    m.m13 = p_rl.y();
    m.m2 = R_rl(2, 0);
    m.m6 = R_rl(2, 1);
    m.m10 = R_rl(2, 2);
    m.m14 = p_rl.z();

    return m;
}

void render_craft(const Eigen::MatrixXd);

// wraps the intial window sizing and framerate
void init_rendering();

// initializes the global camera
void init_camera(const Eigen::Vector3d &pos, const Eigen::Vector3d &target,
                 const Eigen::Vector3d &up, double fovy);

// begin and end for rendering, wraps BeginMode, clearing background etc.
// end_3d finishes the 3d part so 2d things (like the hud) can be drawn on top before end
// (ground draws the grid, which doesn't make sense in space)
void begin(Color background = RAYWHITE, bool ground = true);
void end_3d();
void end();

// text overlay of all the sensor readings next to the original noise free value
void draw_sensor_hud(const Body &b);

} // namespace BSIM::rendering
