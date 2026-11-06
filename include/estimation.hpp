#pragma once
#include "body_state.hpp"
#include <Eigen/Dense>
#include <optional>
#include <type_traits>

namespace BSIM {

// --- state estimation ---
// the two continuous-discrete estimators from chapter 8 (algorithms 1 and 2)
// calling will give one T_out sample (so lines 3 onwards for both)
//
// x_hat (and P for the ekf) are updated each call
// and T_out should be faster than the sensors (line 2)
//
// NX = number of states, NY = size of the measurement from sensor i
// N  = number of euler steps the model takes per T_out
// y  = the new measurement from sensor i, it's std::nullopt if there wasn't a new reading
//
// f, h and their jacobians need to return matrices by value

template <int R, int C> using EstMatrix = Eigen::Matrix<double, R, C>;

// algorithm 1, continuous-discrete observer
// f(x_hat, u) -> x_hat_dot,
// h(x_hat) returns what sensor i should be reading, L is the observer gain for the ith sensor
template <int NX, int NY, typename U, typename F, typename H>
void continuous_discrete_observer(StateVector<NX> &x_hat, const U &u, F &&f, double T_out, int N,
                                  const std::optional<std::type_identity_t<StateVector<NY>>> &y,
                                  H &&h, const EstMatrix<NX, NY> &L) {
    // integrates the the state equation (lines 4-6)
    for (int i = 0; i < N; ++i) {
        x_hat = x_hat + (T_out / N) * f(x_hat, u);
    }

    // measurement update (lines 7-9)
    if (y) {
        x_hat = x_hat + L * (*y - h(x_hat));
    }
}

// algorithm 2, continuous-discrete extended kalman filter
// f(x_hat, u) gives x_hat_dot, and df_dx(x_hat, u) gives A (NX x NX)
// h(x_hat, u) gives what the ith sensor should be reading, dh_dx(x_hat, u) gives C_i (NY x NX)
// Q is the process noise covariance, and R the ith sensors' noise covariance
template <int NX, int NY, typename U, typename F, typename DF, typename H, typename DH>
void continuous_discrete_ekf(StateVector<NX> &x_hat, EstMatrix<NX, NX> &P, const U &u, F &&f,
                             DF &&df_dx, const EstMatrix<NX, NX> &Q, double T_out, int N,
                             const std::optional<std::type_identity_t<StateVector<NY>>> &y, H &&h,
                             DH &&dh_dx, const EstMatrix<NY, NY> &R) {
    // the prediction step (lines 4-8)
    for (int i = 0; i < N; ++i) {
        x_hat = x_hat + (T_out / N) * f(x_hat, u);
        const EstMatrix<NX, NX> A = df_dx(x_hat, u);
        P = P + (T_out / N) * (A * P + P * A.transpose() + Q);
    }

    // measurement update (lines 9-14)
    if (y) {
        const EstMatrix<NY, NX> C = dh_dx(x_hat, u);
        const EstMatrix<NX, NY> L = P * C.transpose() * (R + C * P * C.transpose()).inverse();
        P = (EstMatrix<NX, NX>::Identity(P.rows(), P.cols()) - L * C) * P;
        x_hat = x_hat + L * (*y - h(x_hat, u));
    }
}

} // namespace BSIM
