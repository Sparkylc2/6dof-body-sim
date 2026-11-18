#include "atmosphere.hpp"
#include <algorithm>
#include <cmath>

namespace BSIM {

namespace {

// one layer starting at H_base, with a lapse rate L (dT/dH, negative means colder going up)
struct Layer {
    double H_base, L;
};

constexpr Layer LAYERS[] = {
    {0.0, -0.0065},
    {11000.0, 0.0},
    {20000.0, 0.001},
};
constexpr double H_TOP = 32000.0;

// hydrostatic dp/dH = -rho g0 with p = rho R T, integrated dH up from the base of a layer that
// starts at (T_base, p_base), gives the temperature and pressure there
void through_layer(double L, double dH, double T_base, double p_base, double &T, double &p) {
    T = T_base + L * dH;

    if (L == 0.0) {
        // constant temperature
        p = p_base * std::exp(-ISA_G0 * dH / (ISA_R * T_base));
    } else {
        // linear temperature
        p = p_base * std::pow(T / T_base, -ISA_G0 / (L * ISA_R));
    }
}

} // namespace

double geopotential_altitude(double h) { return ISA_R_EARTH * h / (ISA_R_EARTH + h); }

Atmosphere isa_geopotential(double H) {
    H = std::min(H, H_TOP);

    // climbs up through the layers below H, so each layer's base comes from the one under it
    // (and the profile is continuous by construction, instead of using rounded table values)
    double T = ISA_T0, p = ISA_P0;
    const int n = static_cast<int>(std::size(LAYERS));

    for (int i = 0; i < n; i++) {
        const bool last = i == n - 1 || H < LAYERS[i + 1].H_base;
        const double top = last ? H : LAYERS[i + 1].H_base;
        through_layer(LAYERS[i].L, top - LAYERS[i].H_base, T, p, T, p);

        if (last) {
            break;
        }
    }

    Atmosphere atm;
    atm.T = T;
    atm.p = p;
    atm.rho = p / (ISA_R * T);
    atm.a = std::sqrt(ISA_GAMMA * ISA_R * T);
    return atm;
}

Atmosphere isa(double h) { return isa_geopotential(geopotential_altitude(h)); }

} // namespace BSIM
