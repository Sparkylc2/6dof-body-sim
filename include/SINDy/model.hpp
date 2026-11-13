#pragma once

#include "SINDy/library.hpp"
#include "SINDy/stlsq.hpp"
#include <Eigen/Dense>

namespace BSIM::SINDy {
using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using Eigen::Vector4d;
using Eigen::VectorXd;

struct SindyModel {
    PolynomialLibrary lib;
    MatrixXd Xi; // l x n, in physical units

    void fit(const MatrixXd &X, const MatrixXd &Xdot, double lam) {
        Xi = stlsq(lib.build(X), Xdot, lam);
    }

    VectorXd f(const VectorXd &x) const { return (lib.eval(x) * Xi).transpose(); }

    void print() const {
        const auto names = lib.names();
        for (int k = 0; k < Xi.cols(); k++) {

            std::printf("d(%s)/dt =", lib.state_names[k].c_str());

            for (int j = 0; j < Xi.rows(); j++) {
                if (Xi(j, k) != 0.0) {
                    std::printf(" %+.5g %s", Xi(j, k), names[j].c_str());
                }
            }

            std::printf("\n");
        }
    }
};
} // namespace BSIM::SINDy
