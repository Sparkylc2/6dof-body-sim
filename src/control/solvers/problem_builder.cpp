#include "control/solvers/problem_builder.hpp"

namespace BSIM::Solvers {

Var ProblemBuilder::add_var(int size) {
    // starts wherever the last one ended
    const Var v{n_, size};
    n_ += size;
    return v;
}

void ProblemBuilder::add_quadratic(Var v, const MatrixXd &Q, const VectorXd &target) {
    // expanding it out gives
    // 1/2 (v - t)^T Q (v - t) = 1/2 v^T Q v - (Q t)^T v + const
    // so it's a Q block in P, and a -Q t in q (the const doesn't change the answer so it's dropped)
    P_.push_back({Q, v});
    add_linear(v, -Q * target);
}

void ProblemBuilder::add_linear(Var v, const VectorXd &c) { q_.push_back({c, v}); }

// these just store things, all the actual building happens in solve()
void ProblemBuilder::add_eq(const std::vector<Term> &lhs, const VectorXd &rhs) {
    eq_.push_back({lhs, rhs});
}

void ProblemBuilder::add_le(const std::vector<Term> &lhs, const VectorXd &rhs) {
    le_.push_back({lhs, rhs});
}

void ProblemBuilder::add_norm_le(Var v, double limit) { norms_.push_back({v, limit}); }

ConicResult ProblemBuilder::solve(const ConicSettings &s) const {
    // --- cost ---
    // dense for now, gets made sparse right at the end
    MatrixXd P = MatrixXd::Zero(n_, n_);
    VectorXd q = VectorXd::Zero(n_);

    // each block goes on the diagonal at its variable's spot
    for (const Term &t : P_) {
        P.block(t.v.start, t.v.start, t.v.size, t.v.size) += t.M;
    }
    for (const Term &t : q_) {
        q.segment(t.v.start, t.v.size) += t.M;
    }

    // --- constraints ---
    // counts the rows first so A and b can be sized
    int m = 0;
    for (const auto &c : eq_)
        m += static_cast<int>(c.rhs.size());
    for (const auto &c : le_)
        m += static_cast<int>(c.rhs.size());
    for (const auto &c : norms_)
        m += 1 + c.v.size; // the limit plus the vector

    MatrixXd A = MatrixXd::Zero(m, n_);
    VectorXd b = VectorXd::Zero(m);
    int row = 0; // next free row

    // writes one constraint into A and b
    // each term's M goes in the columns that belong to its variable
    auto put = [&](const Constraint &c) {
        const int rows = static_cast<int>(c.rhs.size());

        for (const Term &t : c.lhs) {
            A.block(row, t.v.start, rows, t.v.size) += t.M;
        }

        b.segment(row, rows) = c.rhs;
        row += rows;
    };

    // clarabel needs the rows in cone order (eq, then le, then norms), so they are stored
    // separately, so they can be added in any order
    ConicProblem prob;

    // zero cone, b - A x = 0
    for (const auto &c : eq_)
        put(c);
    prob.n_eq = row;

    // nonnegative cone, b - A x >= 0 (so A x <= b)
    for (const auto &c : le_)
        put(c);
    prob.n_ineq = row - prob.n_eq;

    // second order cones, b - A x = (limit, v) and the cone says |v| <= limit
    // so the first row is js the limit in b, and the rest are -I on v (since b - (-I v) = v)
    for (const auto &c : norms_) {
        b(row) = c.limit;
        A.block(row + 1, c.v.start, c.v.size, c.v.size) = -MatrixXd::Identity(c.v.size, c.v.size);
        row += 1 + c.v.size;
        prob.soc.push_back(1 + c.v.size);
    }

    // clarabel only reads the upper triangle of P, so giving it the full thing is fine
    prob.P = P.sparseView();
    prob.q = q;
    prob.A = A.sparseView();
    prob.b = b;
    return solve_conic(prob, s);
}

} // namespace BSIM::Solvers
