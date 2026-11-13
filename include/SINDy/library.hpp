#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

namespace BSIM::SINDy {
using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using Eigen::VectorXd;

struct PolynomialLibrary {
    // for now {"du", "dw", "q", "dth"}
    std::vector<std::string> state_names;
    int size() const { return 1 + n() + n() * (n() + 1) / 2; }
    int n() const { return state_names.size(); }

    // takes a sample and returns [1.0, x1, x2, ... x1*x1, x1*x2, ...]
    RowVectorXd eval(const VectorXd &x) const {

        RowVectorXd row(size());

        // write cursor
        int k = 0;

        // constant terms [1.0]
        row(k++) = 1.0;

        // linear terms [x1, x2, x3, ...]
        for (int i = 0; i < n(); i++) {
            row(k++) = x(i);
        }

        // quadratic terms [x1*x1, x1*x2, x1*x3, ...]
        for (int i = 0; i < n(); i++) {
            for (int j = i; j < n(); j++) {
                row(k++) = x(i) * x(j);
            }
        }
        return row;
    };

    // builds Theta(X) (so m x n -> m x size())
    MatrixXd build(const MatrixXd &X) const {

        MatrixXd Theta(X.rows(), size());

        // does [1.0, x1, x2, ..., x1*x1, x1*x2, ...]
        // for each measrement
        for (int m = 0; m < X.rows(); m++) {
            Theta.row(m) = eval(X.row(m).transpose());
        }

        return Theta;
    }

    std::vector<std::string> names() const {
        // constant
        std::vector<std::string> out{"1"};

        // linear
        for (const auto &s : state_names) {
            out.push_back(s);
        }

        // quadratic
        for (int i = 0; i < n(); i++) {
            for (int j = i; j < n(); j++) {
                out.push_back(state_names[i] + "*" + state_names[j]);
            }
        }

        return out;
    }
};
} // namespace BSIM::SINDy
