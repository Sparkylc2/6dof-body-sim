#pragma once
#include <Eigen/Dense>
#include <random>
namespace BSIM {
using Eigen::Vector2d;
using Eigen::Vector3d;

// air density the book uses for all its aircraft (appendix E)
constexpr double RHO_DEFAULT = 1.2682;

// --- general atmosphere ---
// what the air looks like to a body during one step
struct AirState {
    double rho = RHO_DEFAULT;
    Vector3d wind_i = Vector3d::Zero(); // the steady wind, in inertial frame
    Vector3d gust_b = Vector3d::Zero(); // the turbulent wind, in body frame
};

// gusts only need to return the resulting wind velocity through step (in body frame)
struct Gust {
    // could pre-emptively add time as an input but lets see
    virtual Vector3d step(double Va, double dt, std::mt19937 &rng) = 0;
    virtual ~Gust() = default;
};

// the shared atmosphere environment
struct AirSystem {
    double rho = RHO_DEFAULT;
    Vector3d wind_i = Vector3d::Zero(); // shared steady wind
    std::mt19937 rng{0};                // fixed for now
    bool turbulence = true;
};

// --- Dryden atmosphere model ---
// Dryden parameters from Table 4.1
// defaults to low altitude (50 m), light turbulence
struct DrydenParams {
    double sigma_u = 1.06, sigma_v = 1.06, sigma_w = 0.7;
    double L_u = 200.0, L_v = 200.0, L_w = 50.0;
};

// an aircraft's gust filters, but state only, no rng
struct DrydenGust : Gust {
    DrydenParams params;

    double xu = 0.0;   // first order filter state
    Vector2d xv{0, 0}; // second order filter states
    Vector2d xw{0, 0}; // (2 components since state is x, x_dot)

    // advances all three filters in one step and returns
    // [u_g, v_g, w_g] in the body frame
    // the dryden transfer functions from chapter 4, driven by unit white noise
    Vector3d step(double Va, double dt, std::mt19937 &rng) override;
};

} // namespace BSIM
