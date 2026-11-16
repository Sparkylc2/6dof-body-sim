#pragma once
#include "control/solvers/conic.hpp"
#include "control/solvers/problem_builder.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace BSIM::SCvx {
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Solvers::ProblemBuilder;
using Solvers::Var;

// defines the settings (crazy ik)
struct Settings {
    int N = 7;        // num nodes
    double t_h = 6.0; // horizon (dt = t_h / (N-1), [eqn. 7])
    double dt() const { return t_h / (N - 1); }

    double w_ep = 1e3;   // virtual control penalty [eqn. 13]
    double w_prox = 1.0; // proximal weight [eqn. 15]
    double eps = 1e-4;   // allowed y growth per sub-interval [eqn. 11]

    int j_max = 40;         // max scvx iterations
    double tol_prox = 1e-9; // stops once J_prox is below this
    int substeps = 50;      // rk4 steps per sub-interval

    // U [eq 12d]. box is empty for no box, INF norm if theres no norm limit
    VectorXd u_min, u_max;
    double u_norm_max = std::numeric_limits<double>::infinity();

    // convex terminal pieces (cost 1/2 |xi_N - xi_h|^2_Qh, and xi_N[i] = xi_h[i] for i in fix)
    MatrixXd Qh;
    VectorXd xi_h;
    std::vector<int> fix_terminal;

    Solvers::ConicSettings solver;
};

// state and control at each node (x and u)
struct Trajectory {
    std::vector<VectorXd> x, u;
};

// derivatives of smthn w respect to the state and control
struct Jacobians {
    MatrixXd dx; // rows = outputs, cols = state
    MatrixXd du; // rows = outputs, cols = controls
};

// linear model of a sub-interval [eqn. 16b]
// x_{k+1} = A x_k + Bm u_k + Bp u_{k+1} + w
struct Interval {
    MatrixXd A, Bm, Bp;
    VectorXd w;
};

// the variables of [eqn. 16]
struct Vars {
    std::vector<Var> x, u;       // every node
    std::vector<Var> mu_p, mu_m; // virtual control, every sub-interval
};

// what happened in one scvx iteration
struct Iteration {
    double defect = 0;                   // sum |x_{k+1} - f_k|_1 of the point linearised about
    double J_l = 0, J_y = 0, J_prox = 0; // convergence metrics [eqn. 17]
    double l_N = 0;                      // the cost the subproblem found
    Solvers::ConicStatus status = Solvers::ConicStatus::Failed;
};

// the final plan, plus every iteration that got there
struct Result {
    Trajectory Z;
    double defect = 0; // of the final plan (0 means it obeys the real dynamics)
    std::vector<Iteration> history;
};

// definition of the problem [eqn. 1-3]
struct OCP {
    int n_xi = 0; // physical state size
    int n_u = 0;  // control size

    // xi' = F(t, xi, u) [eqn. 1c]
    virtual VectorXd F(double t, const VectorXd &xi, const VectorXd &u) const = 0;
    virtual Jacobians F_jac(double t, const VectorXd &xi, const VectorXd &u) const = 0;

    // running cost L(t, xi, u) [eqn. 3a]
    virtual double L(double t, const VectorXd &xi, const VectorXd &u) const = 0;
    virtual Jacobians L_jac(double t, const VectorXd &xi, const VectorXd &u) const = 0;

    // path inequality constraints g <= 0 [eqn. 2a]
    // (none unless overridden cuz not all problems require constraints)
    virtual VectorXd g(double /*t*/, const VectorXd & /*xi*/, const VectorXd & /*u*/) const {
        return VectorXd(0);
    }
    virtual Jacobians g_jac(double /*t*/, const VectorXd & /*xi*/, const VectorXd & /*u*/) const {
        return {MatrixXd(0, n_xi), MatrixXd(0, n_u)};
    }

    // path equality constraints h == 0 [eqn. 2b]
    virtual VectorXd h(double /*t*/, const VectorXd & /*xi*/, const VectorXd & /*u*/) const {
        return VectorXd(0);
    }
    virtual Jacobians h_jac(double /*t*/, const VectorXd & /*xi*/, const VectorXd & /*u*/) const {
        return {MatrixXd(0, n_xi), MatrixXd(0, n_u)};
    }

    virtual ~OCP() = default;
};

// the mayer form moves the running cost (the integral in eqn. 3a) into the state as l,
// and the path constraints g <= 0 and h == 0 (eqn. 3c, 3d) into the state as y,
// which we can later integrate alongside everything else
//
// tldr is we gain two new scalar states, l and y, which are easily expressed as derivatives
// ergo

// x = (xi, l, y)
inline int nx(const OCP &p) { return p.n_xi + 2; }
// indices of l and y
inline int il(const OCP &p) { return p.n_xi; }
inline int iy(const OCP &p) { return p.n_xi + 1; }

// the full f function (which is used first in eqn. 9)
inline VectorXd f(const OCP &p, double t, const VectorXd &x, const VectorXd &u) {

    const VectorXd xi = x.head(p.n_xi);

    VectorXd xdot(nx(p));

    // xi' = F [eqn. 4]
    xdot.head(p.n_xi) = p.F(t, xi, u);

    // l' = L [eqn 4]
    xdot(il(p)) = p.L(t, xi, u);

    // y' = 1^T max(g(t, xi(t), u(t)), 0)^2 + 1^T h(t, xi(t), u(t))^2 (elementwise) [eqn. 4]
    //
    // l is a scalar anyway (L returns one number). y could have been a vector (one per constraint),
    // however due to the formulation, if any term of g violates <= 0 (bc of the max), or any term
    // of h violates being equal to 0, then the sum is nonzero, so one y is enough
    xdot(iy(p)) = p.g(t, xi, u).cwiseMax(0.0).squaredNorm() + p.h(t, xi, u).squaredNorm();

    return xdot;
}

inline Jacobians f_jac(const OCP &p, double t, const VectorXd &x, const VectorXd &u) {
    const VectorXd xi = x.head(p.n_xi);
    const Jacobians F = p.F_jac(t, xi, u);
    const Jacobians L = p.L_jac(t, xi, u);
    const Jacobians G = p.g_jac(t, xi, u);
    const Jacobians H = p.h_jac(t, xi, u);

    const VectorXd g_plus = p.g(t, xi, u).cwiseMax(0.0);
    const VectorXd hv = p.h(t, xi, u);

    // rows are [xi' | l' | y'],
    // columns are [xi | l | y] (nothing depends on l or y, so 0)
    Jacobians J;
    J.dx = MatrixXd::Zero(nx(p), nx(p));
    J.du = MatrixXd::Zero(nx(p), p.n_u);

    // dF/d(xi) (being n_xi x n_xi)
    J.dx.topLeftCorner(p.n_xi, p.n_xi) = F.dx;
    // dF/du (being n_xi x n_u)
    J.du.topRows(p.n_xi) = F.du;

    // dL/d(xi) (being 1 x n_xi)
    J.dx.row(il(p)).head(p.n_xi) = L.dx;
    // dL/du (being 1 x n_u)
    J.du.row(il(p)) = L.du;

    // d/dz (1' max(g, 0).^2 + 1' h.^2) = 2 max(g, 0)' dg/dz + 2 h' dh/dz
    // (an empty g or h js gives zeroes)
    J.dx.row(iy(p)).head(p.n_xi) = 2.0 * g_plus.transpose() * G.dx + 2.0 * hv.transpose() * H.dx;
    J.du.row(iy(p)) = 2.0 * g_plus.transpose() * G.du + 2.0 * hv.transpose() * H.du;

    return J;
}

// everything integrated along a sub-interval
struct ShotState {
    VectorXd x;      // state [eqn. 9]
    MatrixXd Phi_x;  // dx(t)/dx_k
    MatrixXd Phi_um; // dx(t)/du_k
    MatrixXd Phi_up; // dx(t)/du_{k+1}
};

struct Shot {
    VectorXd x_end; // f_k(x_k, u_k, u_{k+1}) [eqn. 10]
    Interval lin;   // A_k, B_k^-, B_k^+, w_k (the linearization)
};

// s + h * d, for each part of the current shot state
inline ShotState advance(const ShotState &s, const ShotState &d, double h) {
    return {s.x + h * d.x, s.Phi_x + h * d.Phi_x, s.Phi_um + h * d.Phi_um, s.Phi_up + h * d.Phi_up};
}

inline Shot shoot(const OCP &p, double tk, double dt, const VectorXd &xk, const VectorXd &uk,
                  const VectorXd &uk1, int substeps) {

    // the rates of everything in ShotState
    auto rate = [&](double t, const ShotState &s) {
        const double lambda_p = (t - tk) / dt; // FOH weights [eqn. 8]
        const double lambda_m = 1.0 - lambda_p;

        const VectorXd u = lambda_m * uk + lambda_p * uk1; // FOH
        const Jacobians J = f_jac(p, t, s.x, u);           // A(t), B(t) along curr shot

        ShotState d;
        d.x = f(p, t, s.x, u);                        // [eqn. 9]
        d.Phi_x = J.dx * s.Phi_x;                     // Phi_x'  = A Phi_x
        d.Phi_um = J.dx * s.Phi_um + J.du * lambda_m; // Phi_u-' = A Phi_u- + B lam-
        d.Phi_up = J.dx * s.Phi_up + J.du * lambda_p; // Phi_u+' = A Phi_u+ + B lam+
        return d;
    };

    // initial cond are start at x_k, Phi_x = I, and Phi_u = 0
    const int n = nx(p), m = p.n_u;
    ShotState s{xk, MatrixXd::Identity(n, n), MatrixXd::Zero(n, m), MatrixXd::Zero(n, m)};

    // rk4
    const double h = dt / substeps;
    for (int i = 0; i < substeps; i++) {
        const double t = tk + i * h;
        const ShotState k1 = rate(t, s);
        const ShotState k2 = rate(t + h / 2, advance(s, k1, h / 2));
        const ShotState k3 = rate(t + h / 2, advance(s, k2, h / 2));
        const ShotState k4 = rate(t + h, advance(s, k3, h));
        s = advance(advance(advance(advance(s, k1, h / 6), k2, h / 3), k3, h / 3), k4, h / 6);
    }

    Shot out;
    out.x_end = s.x;
    out.lin.A = s.Phi_x;
    out.lin.Bm = s.Phi_um;
    out.lin.Bp = s.Phi_up;
    // makes the linear model exact at (x_k, u_k, u_{k+1})
    out.lin.w = s.x - s.Phi_x * xk - s.Phi_um * uk - s.Phi_up * uk1;
    return out;
}

// shoots every sub-interval along Z
inline std::vector<Interval> linearise(const OCP &p, const Settings &s, double tc,
                                       const Trajectory &Z, double &defect) {
    std::vector<Interval> model;
    defect = 0.0;
    for (int k = 0; k < s.N - 1; k++) {
        const Shot sh = shoot(p, tc + k * s.dt(), s.dt(), Z.x[k], Z.u[k], Z.u[k + 1], s.substeps);
        model.push_back(sh.lin);
        defect += (Z.x[k + 1] - sh.x_end).lpNorm<1>();
    }
    return model;
}

// a 1 x n row with a 1 at i (picks one entry out of a vector, like E_y)
inline MatrixXd pick(int n, int i) {
    MatrixXd E = MatrixXd::Zero(1, n);
    E(0, i) = 1.0;
    return E;
}

inline Trajectory solve_subproblem(const OCP &p, const Settings &s, const VectorXd &xi_c,
                                   const std::vector<Interval> &model, const Trajectory &Zj,
                                   Iteration &it) {
    const int N = s.N, n = nx(p), m = p.n_u;
    const MatrixXd I = MatrixXd::Identity(n, n);
    const MatrixXd Iu = MatrixXd::Identity(m, m);

    ProblemBuilder pb;
    Vars v;
    for (int k = 0; k < N; k++) {
        v.x.push_back(pb.add_var(n));
        v.u.push_back(pb.add_var(m));
    }
    for (int k = 0; k < N - 1; k++) {
        v.mu_p.push_back(pb.add_var(n));
        v.mu_m.push_back(pb.add_var(n));
    }

    // [eqn. 16a] cost
    // l_N (the integrated running cost, linear so exact)
    pb.add_linear(v.x[N - 1], pick(n, il(p)).transpose());
    // convex terminal cost on xi_N
    if (s.Qh.size() != 0) {
        MatrixXd Q = MatrixXd::Zero(n, n);
        Q.topLeftCorner(p.n_xi, p.n_xi) = s.Qh;
        VectorXd target = VectorXd::Zero(n);
        target.head(p.n_xi) = s.xi_h;
        pb.add_quadratic(v.x[N - 1], Q, target);
    }
    // w_ep 1'(mu+ + mu-)
    for (int k = 0; k < N - 1; k++) {
        pb.add_linear(v.mu_p[k], VectorXd::Constant(n, s.w_ep));
        pb.add_linear(v.mu_m[k], VectorXd::Constant(n, s.w_ep));
    }
    // w_prox/2 (|x_k - x_k^j|^2 + |u_k - u_k^j|^2)
    for (int k = 0; k < N; k++) {
        pb.add_quadratic(v.x[k], s.w_prox * I, Zj.x[k]);
        pb.add_quadratic(v.u[k], s.w_prox * Iu, Zj.u[k]);
    }

    // --- linearised dynamics w virtual control [eqn. 16b]----
    //   A x_k + B- u_k + B+ u_{k+1} + mu+ - mu- - x_{k+1} = -w
    for (int k = 0; k < N - 1; k++) {
        const Interval &d = model[k];
        pb.add_eq({{d.A, v.x[k]},
                   {d.Bm, v.u[k]},
                   {d.Bp, v.u[k + 1]},
                   {I, v.mu_p[k]},
                   {-I, v.mu_m[k]},
                   {-I, v.x[k + 1]}},
                  -d.w);
    }

    // ---- y can't grow by more than eps over a sub-interval [eqn. 16c]----
    for (int k = 0; k < N - 1; k++) {
        const MatrixXd Ey = pick(n, iy(p));
        pb.add_le({{Ey, v.x[k + 1]}, {-Ey, v.x[k]}}, VectorXd::Constant(1, s.eps));
    }

    // ----  mu >= 0, written as -mu <= 0 [eqn. 16d] ----
    for (int k = 0; k < N - 1; k++) {
        pb.add_le({{-I, v.mu_p[k]}}, VectorXd::Zero(n));
        pb.add_le({{-I, v.mu_m[k]}}, VectorXd::Zero(n));
    }

    // ---- u_k in U [eqn. 16e]----
    for (int k = 0; k < N; k++) {
        if (s.u_max.size() != 0)
            pb.add_le({{Iu, v.u[k]}}, s.u_max); //  u <= u_max
        if (s.u_min.size() != 0)
            pb.add_le({{-Iu, v.u[k]}}, -s.u_min); // -u <= -u_min
        if (std::isfinite(s.u_norm_max))
            pb.add_norm_le(v.u[k], s.u_norm_max);
    }

    // ---- starts at x_c = (xi_c, l = 0, y = 0) [eqn. 16f]----
    VectorXd xc = VectorXd::Zero(n);
    xc.head(p.n_xi) = xi_c;
    pb.add_eq({{I, v.x[0]}}, xc);

    // ---- terminal equalities [eqn 12f] ----
    for (int i : s.fix_terminal) {
        pb.add_eq({{pick(n, i), v.x[N - 1]}}, VectorXd::Constant(1, s.xi_h(i)));
    }

    // solves and reads everything back
    const Solvers::ConicResult r = pb.solve(s.solver);
    it.status = r.status;
    if (!r.ok())
        return Zj;

    Trajectory Z = Zj;
    for (int k = 0; k < N; k++) {
        Z.x[k] = pb.value(r, v.x[k]);
        Z.u[k] = pb.value(r, v.u[k]);
        it.J_prox += (Z.x[k] - Zj.x[k]).squaredNorm() + (Z.u[k] - Zj.u[k]).squaredNorm(); // 17c
    }

    for (int k = 0; k < N - 1; k++) {
        const VectorXd mu = pb.value(r, v.mu_p[k]) + pb.value(r, v.mu_m[k]);
        it.J_l += mu(il(p)); // [eqn. 17a]
        it.J_y += mu(iy(p)); // [eqn. 17b]
    }
    it.l_N = Z.x[N - 1](il(p));
    return Z;
}

// the scvx loop [eqn. 15]
// linearises around the current plan Z^j, solves the convex subproblem for Z^{j+1},
// then repeats until J_prox is small or j_max runs out
inline Result solve(const OCP &p, const Settings &s, double tc, const VectorXd &xi_c,
                    Trajectory Z) {
    Result res;
    for (int j = 0; j < s.j_max; j++) {
        Iteration it;
        const std::vector<Interval> model = linearise(p, s, tc, Z, it.defect); // around Z^j
        Z = solve_subproblem(p, s, xi_c, model, Z, it);                        // Z^{j+1}
        res.history.push_back(it);

        const bool solved = it.status == Solvers::ConicStatus::Solved ||
                            it.status == Solvers::ConicStatus::AlmostSolved;
        if (!solved || it.J_prox < s.tol_prox) {
            break;
        }
    }

    // final shoot once more to see how far the final plan is from the real dynamics
    linearise(p, s, tc, Z, res.defect);
    res.Z = Z;
    return res;
}

// ---- nmpc [section 3, fig. 1] ----

// the guess for the next horizon (shifts the last plan along by K nodes)
// the K new nodes on the end hold u_N and integrate on from the node before
inline Trajectory warm_start(const OCP &p, const Settings &s, double tc, const Trajectory &Z,
                             int K) {
    const int N = s.N;
    Trajectory W = Z;

    // shift along by K nodes
    for (int k = 0; k < N - K; k++) {
        W.x[k] = Z.x[k + K];
        W.u[k] = Z.u[k + K];
    }

    // the K new nodes, u is kept at u_N
    for (int k = N - K; k < N; k++) {
        W.u[k] = Z.u[N - 1];
        const double t0 = tc + s.t_h + (k - (N - K)) * s.dt();
        VectorXd x = W.x[k - 1];
        const double h = s.dt() / s.substeps;

        for (int i = 0; i < s.substeps; i++) {
            const double t = t0 + i * h;
            const VectorXd k1 = f(p, t, x, W.u[k]);
            const VectorXd k2 = f(p, t + h / 2, x + h / 2 * k1, W.u[k]);
            const VectorXd k3 = f(p, t + h / 2, x + h / 2 * k2, W.u[k]);
            const VectorXd k4 = f(p, t + h, x + h * k3, W.u[k]);
            x += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
        }
        W.x[k] = x;
    }

    return W;
}

// the control the plan gives at time t, foh between the nodes [eqn. 8]
// (clamped to the horizon, so past the end it holds u_N)
inline VectorXd control_at(const Trajectory &Z, double tc, double dt, double t) {
    const int N = static_cast<int>(Z.u.size());
    const double s = std::clamp((t - tc) / dt, 0.0, static_cast<double>(N - 1));
    const int k = std::min(static_cast<int>(s), N - 2);

    const double lambda_p = s - k;

    return (1.0 - lambda_p) * Z.u[k] + lambda_p * Z.u[k + 1];
}

} // namespace BSIM::SCvx
