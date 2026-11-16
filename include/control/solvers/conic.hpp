#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

namespace BSIM::Solvers {
using Eigen::VectorXd;
using SparseMatrixd = Eigen::SparseMatrix<double, Eigen::ColMajor>;

// --- convex conic program (clarabel's standard form) ---
//
// (icl as much as I wanted to write the solver myself it's out of scope for this project)
// this structure is AI generated, however I've read through it (as you can probably tell)
//
// takes the form (standard quadratic):
//      min  1/2 x^T P x + q^T x
//      s.t. Ax + s = b, s in K (so b - Ax = s)
//
// s is a slack, so each block of rows says "b - A z is in this cone"
//      Ax = b      (an equality constraint, s = 0)
//      Ax <= b     (an inequality, s >= 0)
//      ||v|| <= t  (norm limit, s = (t, v))

// K is made of these blocks, and the rows of A and b have to be in the order:
//      zero cone (n_eq rows):             Ax = b (equality constraints, s = 0)
//      nonnegative cone (n_ineq rows):    Ax <= b (inequality constraints (s >= 0)
//      second order cones (soc[i] rows each):
//                                          s = (t, v), |v| <= t (so norm limits like |u| <= u_max)
struct ConicProblem {
    SparseMatrixd P; // n x n, only the upper triangle is read
    VectorXd q;      // n
    SparseMatrixd A; // m x n
    VectorXd b;      // m
    int n_eq = 0;
    int n_ineq = 0;
    std::vector<int> soc; // size of each second order cone (the t plus the vector)
};

struct ConicSettings {
    int max_iter = 200;
    double tol = 1e-8; // used for the gap (abs and rel) and the feasibility tolerances
    bool verbose = false;
};

enum class ConicStatus {
    Solved,
    AlmostSolved,
    PrimalInfeasible, // no x satisfies all the constraints at once
    DualInfeasible,   // the cost is unbounded (can keep going down forever)
    MaxIter,
    Failed, // numerical trouble, no progress etc
};

struct ConicResult {
    VectorXd x; // solution
    VectorXd y; // the lagrange multipliers of the rows (so each constraint) of A x + s = b
    VectorXd s; // slack
    ConicStatus status = ConicStatus::Failed;
    int iters = 0;
    double solve_time = 0.0;

    bool ok() const { return status == ConicStatus::Solved || status == ConicStatus::AlmostSolved; }
};

// solves it with clarabel (an interior point method)
ConicResult solve_conic(const ConicProblem &prob, const ConicSettings &s = {});

const char *to_string(ConicStatus s);

} // namespace BSIM::Solvers
