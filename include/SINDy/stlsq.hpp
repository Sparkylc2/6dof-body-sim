#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

namespace BSIM::SINDy {
using Eigen::MatrixXd;
using Eigen::Vector4d;
using Eigen::VectorXd;

using MaskXX = Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>;
using MaskX = Eigen::Array<bool, Eigen::Dynamic, 1>;

// solves Ax = y with the normal x = (A^T A)^-1 A^T y
inline MatrixXd least_squares(const MatrixXd &y, const MatrixXd &A) {
    return (A.transpose() * A).inverse() * A.transpose() * y;
}

// columns of A where keep is true (which is like numpy's theta[:, mask])
static MatrixXd select_columns(const MatrixXd &A, const MaskX &keep) {
    MatrixXd out(A.rows(), keep.count());

    int c = 0;

    for (int j = 0; j < A.cols(); j++) {
        if (keep(j)) {
            out.col(c++) = A.col(j);
        }
    }
    return out;
}

// basically straight from the tutorial repo (rewritten)
// theta is (m x l), x_dot is m x n derivatives, and this returns xi (l x n)
inline MatrixXd stlsq(const MatrixXd &theta, const MatrixXd &x_dot, double lam, int max_iter = 20) {
    // initial least squares fit
    MatrixXd xi = least_squares(x_dot, theta);

    // masks the coefficients and then zeros the ones below the lam threshold
    MaskXX mask = xi.array().abs() >= lam;
    xi.array() *= mask.cast<double>();

    // copy of the mask for the convergence check
    MaskXX prev_mask = mask;

    // least squares fit per state
    for (int state = 0; state < x_dot.cols(); ++state) {
        for (int it = 0; it < max_iter; ++it) {
            // nothing left
            if (mask.col(state).count() == 0)
                break;

            // refit is only run on terms that weren't removed
            const Eigen::VectorXd coeffs =
                least_squares(x_dot.col(state), select_columns(theta, mask.col(state)));

            // put the terms back
            int c = 0;
            for (int j = 0; j < xi.rows(); ++j) {
                if (mask(j, state)) {
                    xi(j, state) = coeffs(c++);
                }
            }

            // updates the mask and sparsify again
            mask.col(state) = xi.col(state).array().abs() >= lam;
            xi.array() *= mask.cast<double>();

            // it's converged if no terms were added or removed
            if ((mask == prev_mask).all()) {
                break;
            }
            prev_mask = mask;
        }
    }
    return xi;
}

} // namespace BSIM::SINDy
