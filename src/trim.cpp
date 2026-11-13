#include "trim.hpp"
#include "body_system.hpp"
#include "frames.hpp"

namespace BSIM {
using Eigen::Matrix3d;

AircraftState LongitudinalTrim::state(const Vector3d &pos) const {
    AircraftState x;
    x.pos() = pos;
    // no wind
    x.vel() = Vector3d(Va * std::cos(alpha), 0.0, Va * std::sin(alpha));
    x.quat() = euler_to_quat(0.0, theta(), 0.0);
    x.rates() = Vector3d::Zero();
    return x;
}

LongitudinalTrim trim_longitudinal(const FixedWing &plane, double Va, double gamma, double rho) {
    LongitudinalTrim trim;
    trim.Va = Va;
    trim.gamma = gamma;

    AirState air;
    air.rho = rho;

    const UAVBodySystem sys;

    // z = [alpha, de, dt] -> r = [u_dot, w_dot, q_dot]
    auto residual = [&](const Vector3d &z) {
        LongitudinalTrim t = trim;
        t.alpha = z(0);
        const AircraftState x = t.state();
        const FixedWingControls u{.de = z(1), .da = 0.0, .dr = 0.0, .dt = z(2)};

        const Wrench w = plane.wrench(x, air, u);
        const auto x_dot = sys.derivatives(plane, x, w.F, w.M);
        return Vector3d(x_dot(IDX::VEL), x_dot(IDX::VEL + 2), x_dot(IDX::RATES + 1));
    };

    Vector3d z(0.05, 0.0, 0.5);
    Vector3d r = residual(z);

    constexpr int max_iter = 50;
    constexpr double tol = 1e-10;
    constexpr double h = 1e-7;

    // newton method
    for (int i = 0; i < max_iter && r.norm() > tol; ++i) {
        Matrix3d J;
        for (int j = 0; j < 3; ++j) {
            Vector3d zp = z;
            zp(j) += h;
            J.col(j) = (residual(zp) - r) / h;
        }
        z = z - J.colPivHouseholderQr().solve(r);
        r = residual(z);
    }

    trim.alpha = z(0);
    trim.controls = {.de = z(1), .da = 0.0, .dr = 0.0, .dt = z(2)};
    trim.residual = r.norm();
    trim.converged = trim.residual <= tol;
    return trim;
}

} // namespace BSIM
