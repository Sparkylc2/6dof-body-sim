#include "force_models.hpp"
#include "air_system.hpp"
#include "frames.hpp"
#include <cmath>

namespace uav {

Wrench gravity(double m, const Vector4d &e) {
    Wrench w;

    Eigen::Vector3d g = {0, 0, m * GRAVITY};
    // gravity is along +d in NED, rotation(e) is body -> vehicle so you need to transpose
    Eigen::Matrix3d R_grav = rotation(e).transpose();

    w.F = R_grav * g;

    return w;
}

AirData air_data(const Vector3d &vel_ground_b, const Vector4d &quat, const AirState &air) {
    // collect the wind info
    const Vector3d wind_b = rotation(quat).transpose() * air.wind_i + air.gust_b;
    const Vector3d vel_air_b = vel_ground_b - wind_b;

    // get the air-relative components
    const double u_r = vel_air_b.x(), v_r = vel_air_b.y(), w_r = vel_air_b.z();

    AirData air_info;

    air_info.Va = vel_air_b.norm();
    air_info.alpha = std::atan2(w_r, u_r);
    // same as asin(v_r / Va) when Va > 0, but this one is still defined (ie 0) at rest
    air_info.beta = std::atan2(v_r, std::sqrt(u_r * u_r + w_r * w_r));
    air_info.rho = air.rho;

    return air_info;
}

} // namespace uav
