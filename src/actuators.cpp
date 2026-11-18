#include "actuators.hpp"
#include <algorithm>
#include <cmath>

namespace BSIM {

double Actuator::step(double cmd, double dt) {
    switch (kind) {
    case Kind::Ideal:
        rate = (std::clamp(cmd, -pos_max, pos_max) - pos) / dt;
        pos = std::clamp(cmd, -pos_max, pos_max);
        return pos;

    case Kind::FirstOrder: {
        // exact for the linear part over each sub step, then goes through the limits
        const int n = std::max(1, static_cast<int>(std::ceil(dt / (0.1 * tau))));
        const double h = dt / n;
        for (int i = 0; i < n; i++) {
            const double target = pos + (cmd - pos) * (1.0 - std::exp(-h / tau));
            rate = std::clamp((target - pos) / h, -rate_max, rate_max);
            pos = std::clamp(pos + rate * h, -pos_max, pos_max);
        }
        return pos;
    }

    case Kind::SecondOrder: {
        // semi implicit euler
        const int n = std::max(1, static_cast<int>(std::ceil(dt * wn / 0.02)));
        const double h = dt / n;

        for (int i = 0; i < n; i++) {
            const double acc = wn * wn * (cmd - pos) - 2.0 * zeta * wn * rate;
            rate = std::clamp(rate + acc * h, -rate_max, rate_max);

            const double next = pos + rate * h;
            pos = std::clamp(next, -pos_max, pos_max);

            if (pos != next) {
                rate = 0.0; // hit the end
            }
        }
        return pos;
    }
    }
    return pos;
}

void Actuator::reset(double p) {
    pos = p;
    rate = 0.0;
}

} // namespace BSIM
