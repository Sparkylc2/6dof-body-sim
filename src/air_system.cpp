#include "air_system.hpp"
#include <algorithm>
#include <cmath>

namespace BSIM {

Vector3d DrydenGust::step(double Va, double dt, std::mt19937 &rng) {
    std::normal_distribution<double> noise(0.0, 1.0);

    // the gains and poles all scale with Va, so i'm clamping it to 1 so they don't collapse onto 0
    Va = std::max(Va, 1.0);
    const DrydenParams &p = params;

    // --- u_g, first order, H_u = sigma_u sqrt(2 Va / L_u) / (s + Va / L_u) ---
    // x_dot = -a x + eta
    const double au = Va / p.L_u;
    xu += dt * (-au * xu) + std::sqrt(dt) * noise(rng);
    const double u_g = p.sigma_u * std::sqrt(2.0 * Va / p.L_u) * xu;

    // --- v, w: second order, H = sigma sqrt(3 Va / L) (s + a / sqrt(3)) / (s + a)^2 ---
    // u can use a standard recast with, x = [x1, x2]
    //   x1_dot = x2
    //   x2_dot = -a^2 x1 - 2a x2 + eta
    //   out    = K (b x1 + x2)
    // stepped with euler-maruyama, the white noise integrated over dt has variance dt so it gets
    // sqrt(dt) not dt (the time constant is L / Va, seconds, so euler should be good enough here)
    auto second_order = [&](Vector2d &x, double sigma, double L) {
        const double a = Va / L;
        const double b = a / std::sqrt(3.0);
        const double K = sigma * std::sqrt(3.0 * Va / L);

        const Vector2d x_dot(x.y(), -a * a * x.x() - 2.0 * a * x.y());
        x += dt * x_dot;
        x.y() += std::sqrt(dt) * noise(rng);
        return K * (b * x.x() + x.y());
    };

    const double v_g = second_order(xv, p.sigma_v, p.L_v);
    const double w_g = second_order(xw, p.sigma_w, p.L_w);

    return {u_g, v_g, w_g};
}

} // namespace BSIM
