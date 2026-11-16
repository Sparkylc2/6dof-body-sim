#pragma once
#include "control/solvers/conic.hpp"
#include <Eigen/Dense>

namespace BSIM::Solvers {
using Eigen::MatrixXd;
using Eigen::VectorXd;

// a block of decision variables (eg x_3 or u_1), just remembers where it sits in the big vector
struct Var {
    int start = 0;
    int size = 0;
};

// one piece of a constraint, M * v
struct Term {
    MatrixXd M;
    Var v;
};

class ProblemBuilder {
public:
    Var add_var(int size);

    // cost
    void add_quadratic(Var v, const MatrixXd &Q, const VectorXd &target); // + 1/2 |v - target|^2_Q
    void add_linear(Var v, const VectorXd &c);                            // + c' v

    // constraints, left side is sum of terms
    void add_eq(const std::vector<Term> &lhs, const VectorXd &rhs); // sum M v == rhs
    void add_le(const std::vector<Term> &lhs, const VectorXd &rhs); // sum M v <= rhs
    void add_norm_le(Var v, double limit);                          // |v| <= limit

    ConicResult solve(const ConicSettings &s = {}) const;
    VectorXd value(const ConicResult &r, Var v) const { return r.x.segment(v.start, v.size); }

private:
    struct Constraint {
        std::vector<Term> lhs;
        VectorXd rhs;
    };
    struct NormLimit {
        Var v;
        double limit;
    };

    int n_ = 0;           // total number of variables so far
    std::vector<Term> P_; // quadratic cost blocks (M is Q, v says where it goes)
    std::vector<Term> q_; // linear cost blocks (M is c as a column)
    std::vector<Constraint> eq_, le_;
    std::vector<NormLimit> norms_;
};

} // namespace BSIM::Solvers
