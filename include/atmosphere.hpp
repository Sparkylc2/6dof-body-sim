#pragma once

// --- international standard atmosphere ---
// (ISA, the same as the 1976 US standard atmosphere up to 32 km)
// the air is treated as a perfect gas with a temperature that changes linearly in each layer
//      0 - 11 km:  troposphere, cools at 6.5 K/km
//      11 - 20 km: tropopause, constant temperature
//      20 - 32 km: lower stratosphere, warms at 1 K/km
// past 32 km it's held at the 32 km values
namespace BSIM {

// sea level
constexpr double ISA_T0 = 288.15;    // temperature [K]
constexpr double ISA_P0 = 101325.0;  // pressure [Pa]
constexpr double ISA_RHO0 = 1.225;   // density [kg/m^3]
constexpr double ISA_A0 = 340.29399; // speed of sound [m/s]

constexpr double ISA_R = 287.05287;       // specific gas constant for air [J/(kg K)]
constexpr double ISA_GAMMA = 1.4;         // ratio of specific heats
constexpr double ISA_G0 = 9.80665;        // the gravity the layers are defined with [m/s^2]
constexpr double ISA_R_EARTH = 6356766.0; // radius used for geopotential altitude [m]

// the air at one altitude
struct Atmosphere {
    double T = ISA_T0;     // temperature [K]
    double p = ISA_P0;     // pressure [Pa]
    double rho = ISA_RHO0; // density [kg/m^3]
    double a = ISA_A0;     // speed of sound [m/s]
};

// the layers are defined in geopotential altitude (as if gravity didn't fall off with height),
// this converts from the actual (geometric) altitude above sea level
// (the difference is small, about 13 m at 9 km)
double geopotential_altitude(double h);

// the atmosphere at a geometric altitude h above sea level [m]
// (below sea level it just carries on with the troposphere)
Atmosphere isa(double h);

// the atmosphere at a geopotential altitude H [m] (what matlab's atmosisa takes)
Atmosphere isa_geopotential(double H);

} // namespace BSIM
