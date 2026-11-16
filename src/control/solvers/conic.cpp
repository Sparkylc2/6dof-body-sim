#include "control/solvers/conic.hpp"
#include <clarabel.hpp> // only ever included here, everything else goes through conic.hpp
#include <stdexcept>

namespace BSIM::Solvers {
namespace {

ConicStatus convert(clarabel::SolverStatus s) {
    using clarabel::SolverStatus;
    switch (s) {
    case SolverStatus::Solved:
        return ConicStatus::Solved;
    case SolverStatus::AlmostSolved:
        return ConicStatus::AlmostSolved;
    case SolverStatus::PrimalInfeasible:
    case SolverStatus::AlmostPrimalInfeasible:
        return ConicStatus::PrimalInfeasible;
    case SolverStatus::DualInfeasible:
    case SolverStatus::AlmostDualInfeasible:
        return ConicStatus::DualInfeasible;
    case SolverStatus::MaxIterations:
        return ConicStatus::MaxIter;
    default:
        return ConicStatus::Failed;
    }
}

} // namespace

ConicResult solve_conic(const ConicProblem &prob, const ConicSettings &s) {
    const Eigen::Index n = prob.q.size();
    const Eigen::Index m = prob.b.size();

    int soc_rows = 0;
    for (int d : prob.soc) {
        soc_rows += d;
    }

    if (prob.P.rows() != n || prob.P.cols() != n || prob.A.rows() != m || prob.A.cols() != n ||
        prob.n_eq + prob.n_ineq + soc_rows != m) {
        throw std::invalid_argument("solve_conic - sizes of P, q, A, b and the cones don't match");
    }

    // the cones in row order
    std::vector<clarabel::SupportedConeT<double>> cones;
    if (prob.n_eq > 0) {
        cones.push_back(clarabel::ZeroConeT<double>(prob.n_eq));
    }
    if (prob.n_ineq > 0) {
        cones.push_back(clarabel::NonnegativeConeT<double>(prob.n_ineq));
    }
    for (int d : prob.soc) {
        cones.push_back(clarabel::SecondOrderConeT<double>(d));
    }

    clarabel::DefaultSettings<double> cs = clarabel::DefaultSettings<double>::default_settings();
    cs.max_iter = static_cast<uint32_t>(s.max_iter);
    cs.verbose = s.verbose;
    cs.tol_gap_abs = s.tol;
    cs.tol_gap_rel = s.tol;
    cs.tol_feas = s.tol;

    // clarabel apparently wants compressed column major matrices with only the upper triangle of P
    SparseMatrixd P = prob.P.triangularView<Eigen::Upper>();
    SparseMatrixd A = prob.A;
    P.makeCompressed();
    A.makeCompressed();
    VectorXd q = prob.q;
    VectorXd b = prob.b;

    clarabel::DefaultSolver<double> solver(P, q, A, b, cones, cs);
    solver.solve();
    const clarabel::DefaultSolution<double> sol = solver.solution();

    // copy results out
    ConicResult r;
    r.x = sol.x;
    r.y = sol.z;
    r.s = sol.s;
    r.status = convert(sol.status);
    r.iters = static_cast<int>(sol.iterations);
    r.solve_time = sol.solve_time;
    return r;
}

const char *to_string(ConicStatus s) {
    switch (s) {
    case ConicStatus::Solved:
        return "solved";
    case ConicStatus::AlmostSolved:
        return "almost solved";
    case ConicStatus::PrimalInfeasible:
        return "primal infeasible";
    case ConicStatus::DualInfeasible:
        return "dual infeasible";
    case ConicStatus::MaxIter:
        return "max iterations";
    case ConicStatus::Failed:
        return "failed";
    }
    return "?";
}

} // namespace BSIM::Solvers
