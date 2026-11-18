// scvx nmpc flying a drone round an L shaped course through randomly placed balls
//
//   scvx plans the next t_h seconds in the horizontal plane with a simple turn rate model
//   (constant speed, chi' = r) and the balls as path constraints g <= 0 (eqn. 2a)
//   the lateral autopilot then flies the course the plan gives (with the plan's turn rate as a
//   roll feedforward), and the longitudinal autopilot holds the height and airspeed (the balls
//   are all centred at the flight altitude, so the plan only has to think in 2d)
//   every dt (one node) it replans from the measured state, warm started from the last plan
//   (algorithm 1, K = 1). in the window the solve runs on a second thread so the frame rate
//   doesn't die, and the plane keeps flying the old plan until the new one is ready
#include "body.hpp"
#include "control/autopilots/autopilot.hpp"
#include "control/autopilots/lateral_autopilot.hpp"
#include "control/autopilots/longitudinal_autopilot.hpp"
#include "control/scvx/scvx.hpp"
#include "frames.hpp"
#include "render/render.hpp"
#include "scenes/scenes.hpp"
#include "sensors.hpp" // wrap_angle, DEG
#include "trim.hpp"
#include "world.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <future>
#include <memory>
#include <numbers>
#include <random>
#include <raylib.h>
#include <rlgl.h>
#include <vector>

namespace BSIM::scenes {
namespace {

using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::VectorXd;

// course
// a polyline of (n, e) corners, past the last corner it carries on straight
struct Course {
    std::vector<Vector2d> pts;

    double length() const {
        double L = 0.0;
        for (size_t i = 0; i + 1 < pts.size(); i++) {
            L += (pts[i + 1] - pts[i]).norm();
        }
        return L;
    }

    // which leg arc length s is on/how far along it
    void leg(double s, size_t &i, double &along) const {
        for (i = 0; i + 2 < pts.size(); i++) {
            const double Li = (pts[i + 1] - pts[i]).norm();
            if (s < Li) {
                break;
            }
            s -= Li;
        }
        along = s;
    }

    Vector2d at(double s) const {
        size_t i;
        double along;
        leg(std::max(s, 0.0), i, along);
        return pts[i] + along * (pts[i + 1] - pts[i]).normalized();
    }

    // course angle of the leg at s (chi = 0 is north, pi/2 is east)
    double heading(double s) const {
        size_t i;
        double along;

        leg(std::max(s, 0.0), i, along);

        const Vector2d d = pts[i + 1] - pts[i];
        return std::atan2(d.y(), d.x());
    }

    // closest point to p on the whole course (the last leg is extended so the carrot always exists)
    // and its arc length

    Vector2d closest(const Vector2d &p, double *s_out = nullptr) const {
        Vector2d best = pts[0];
        double best_d = INF, s_best = 0.0, s0 = 0.0;

        for (size_t i = 0; i + 1 < pts.size(); i++) {
            const Vector2d a = pts[i], d = pts[i + 1] - pts[i];
            const double L = d.norm();
            const bool last = i + 2 == pts.size();
            const double t = std::clamp((p - a).dot(d) / L, 0.0, last ? INF : L);
            const Vector2d c = a + t * d / L;

            if ((p - c).norm() < best_d) {
                best_d = (p - c).norm();
                best = c;
                s_best = s0 + t;
            }
            s0 += L;
        }

        if (s_out) {
            *s_out = s_best;
        }
        return best;
    }
};

// an obstacle, drawn as a ball centred at the flight altitude
struct Ball {
    Vector2d c; // (n, e)
    double R = 0.0;
    int colour = 0;
};

// --- planner model [eqn. 1-3] ---
// xi = (n, e, chi), u = r (the course rate)

// constant ground speed V, so  n' = V cos chi,  e' = V sin chi,  chi' = r
struct TurnRateOCP : SCvx::OCP {
    double V = 18.0;
    double margin = 3.0;  // extra clearance round each ball (covers the tracking error too)
    double w_path = 0.01; // pulls plan back onto course line
    double w_r = 5.0;     // keeps the turns gentle
    const Course *course = nullptr;
    const std::vector<Ball> *balls = nullptr;

    TurnRateOCP() {
        n_xi = 3;
        n_u = 1;
    }

    VectorXd F(double /*t*/, const VectorXd &xi, const VectorXd &u) const override {
        return Vector3d(V * std::cos(xi(2)), V * std::sin(xi(2)), u(0));
    }

    SCvx::Jacobians F_jac(double /*t*/, const VectorXd &xi, const VectorXd & /*u*/) const override {
        SCvx::Jacobians J{MatrixXd::Zero(3, 3), MatrixXd::Zero(3, 1)};
        J.dx(0, 2) = -V * std::sin(xi(2));
        J.dx(1, 2) = V * std::cos(xi(2));
        J.du(2, 0) = 1.0;
        return J;
    }

    // L = w_path |p - closest point on the course|^2 + w_r r^2
    double L(double /*t*/, const VectorXd &xi, const VectorXd &u) const override {
        const Vector2d p = xi.head<2>();
        return w_path * (p - course->closest(p)).squaredNorm() + w_r * u(0) * u(0);
    }

    SCvx::Jacobians L_jac(double /*t*/, const VectorXd &xi, const VectorXd &u) const override {
        const Vector2d p = xi.head<2>();
        SCvx::Jacobians J{MatrixXd::Zero(1, 3), MatrixXd::Zero(1, 1)};

        // the gradient of the squared distance to a line is 2 (p - closest)
        J.dx.leftCols(2) = 2.0 * w_path * (p - course->closest(p)).transpose();
        J.du(0, 0) = 2.0 * w_r * u(0);

        return J;
    }

    // g_i = R_i + margin - |p - c_i| <= 0 (one per ball)
    VectorXd g(double /*t*/, const VectorXd &xi, const VectorXd & /*u*/) const override {
        VectorXd out(balls->size());

        for (size_t i = 0; i < balls->size(); i++) {
            const Ball &b = (*balls)[i];
            out(static_cast<int>(i)) = b.R + margin - (xi.head<2>() - b.c).norm();
        }

        return out;
    }

    SCvx::Jacobians g_jac(double /*t*/, const VectorXd &xi, const VectorXd & /*u*/) const override {
        const int m = static_cast<int>(balls->size());
        SCvx::Jacobians J{MatrixXd::Zero(m, 3), MatrixXd::Zero(m, 1)};

        for (int i = 0; i < m; i++) {
            const Vector2d d = xi.head<2>() - (*balls)[static_cast<size_t>(i)].c;
            J.dx.block(i, 0, 1, 2) = -d.transpose() / std::max(d.norm(), 1e-9);
        }

        return J;
    }
};

// course along the plan at time t
// it's the integral of the foh course rate, so a quadratic between nodes rather than a straight
// line
double plan_course(const SCvx::Trajectory &Z, double t_plan, double dt, double t) {
    const int N = static_cast<int>(Z.x.size());
    const double s = std::clamp((t - t_plan) / dt, 0.0, static_cast<double>(N - 1));
    const int k = std::min(static_cast<int>(s), N - 2);

    const double tau = (s - k) * dt;
    const double uk = Z.u[static_cast<size_t>(k)](0), uk1 = Z.u[static_cast<size_t>(k + 1)](0);

    return Z.x[static_cast<size_t>(k)](2) + uk * tau + (uk1 - uk) * tau * tau / (2.0 * dt);
}

// course rate along the plan at time t (foh, so straight lines between the nodes [eqn. 8])
double plan_rate(const SCvx::Trajectory &Z, double t_plan, double dt, double t) {
    return SCvx::control_at(Z, t_plan, dt, t)(0);
}

// where the plan says the plane should be at time t (integrates the course from the node before)
Vector2d plan_position(const SCvx::Trajectory &Z, double t_plan, double dt, double V, double t) {
    const int N = static_cast<int>(Z.x.size());
    const double s = std::clamp((t - t_plan) / dt, 0.0, static_cast<double>(N - 1));
    const int k = std::min(static_cast<int>(s), N - 2);
    const double tk = t_plan + k * dt;

    const int steps = 8;
    const double h = (t_plan + s * dt - tk) / steps;
    Vector2d p = Z.x[static_cast<size_t>(k)].head<2>();
    for (int i = 0; i < steps; i++) {
        const double chi = plan_course(Z, t_plan, dt, tk + (i + 0.5) * h);
        p += V * h * Vector2d(std::cos(chi), std::sin(chi));
    }
    return p;
}

struct Mission {
    // flight
    double Va = 18.0;  // airspeed
    double h_c = 30.0; // altitude everything happens at
    int n_balls = 8;

    World world;
    FixedWing *body = nullptr;
    LongitudinalTrim trim;
    LateralDesignParams lat;
    LateralAutopilot lat_ap;
    LongitudinalAutopilot lon;

    // planner
    Course course;
    std::vector<Ball> balls;
    TurnRateOCP ocp;
    SCvx::Settings s;
    double lead = 0.3; // how far ahead on the plan the commands are taken from (roll lag) [s]

    // pulls the plane back onto the plan's position between replans
    // chi_c = chi_plan - chi_inf (2 / pi) atan(k_ct e_y), e_y is how far right of the plan it is
    double k_ct = 0.15;
    double chi_inf = 30.0 * DEG;

    SCvx::Trajectory Z; // the current plan
    SCvx::Result last;  // the last solve
    double t_plan = 0.0;
    double solve_ms = 0.0;
    bool have_plan = false;

    // solves on a second thread
    struct Solved {
        SCvx::Result res;
        double ms = 0.0;
    };

    bool async = false;
    std::future<Solved> pending;
    double pending_t = 0.0; // when the pending plan was measured from

    // misc vars
    double t = 0.0;
    double s_progress = 0.0; // arc len along the course
    double closest_pass = INF;
    int hits = 0;
    int run = 0;

    std::vector<Vector3d> trail;
    double next_trail = 0.0;
    std::mt19937 rng{7};

    Mission() {
        auto plane = std::make_unique<FixedWing>(generate_sport_drone());
        body = plane.get();
        world.bodies.bodies.push_back(std::move(plane));
        world.atmosphere.turbulence = false;

        trim = trim_longitudinal(*body, Va);
        body->trim = trim.controls;
        lon.set_trim(trim);
        lon.clock.Ts = body->sensors.params.Ts;

        // same as the plane scene
        lat.da_max = 20.0 * DEG;
        lat.e_phi_max = 80.0 * DEG;
        lat.zeta_phi = 1.4;
        lat.ki_phi = 0.1;
        lat.W_chi = 12.0;
        lat.zeta_chi = 1.0;
        lat.Vg = Va;
        lat.phi_max = 30.0 * DEG;
        lat.dr_max = 30.0 * DEG;
        lat.e_beta_max = 10.0 * DEG;
        lat.zeta_beta = 0.9;
        lat_ap.clock.Ts = body->sensors.params.Ts;
        lat_ap.design(body->params, body->m, body->J, RHO_DEFAULT, Va, lat);

        // north for 550 m, then a right turn and east for 550 m
        course.pts = {Vector2d(-300, -250), Vector2d(250, -250), Vector2d(250, 300)};

        // planner settings, 8 s horizon with a node every second
        s.N = 15;
        s.t_h = 8.0;
        s.substeps = 15;
        s.w_prox = 0.3;

        // turn rate limit from a 25 deg bank, r = g tan(phi) / V (gives the autopilot a lil room)
        const double r_max = 9.81 * std::tan(25.0 * DEG) / Va;
        s.u_min = VectorXd::Constant(1, -r_max);
        s.u_max = VectorXd::Constant(1, r_max);

        // terminal cost pulls the end of the plan to the carrot (V t_h further along the course)
        s.Qh = Vector3d(0.05, 0.05, 20.0).asDiagonal();
        s.xi_h = Vector3d::Zero();

        ocp.V = Va;
        ocp.course = &course;
        ocp.balls = &balls;

        reset();
    }

    // random balls near the course line, spaced so there's always a way through
    void place_balls() {
        balls.clear();
        const double L = course.length();

        std::uniform_real_distribution<double> s_dist(90.0, L - 40.0), off(-12.0, 12.0),
            R_dist(5.0, 10.0);

        int colour = 0;
        for (int tries = 0; static_cast<int>(balls.size()) < n_balls && tries < 2000; tries++) {
            const double sb = s_dist(rng);
            const double chi = course.heading(sb);
            const Vector2d side(-std::sin(chi), std::cos(chi));

            Ball b{course.at(sb) + off(rng) * side, R_dist(rng), colour};
            const bool clear = std::none_of(balls.begin(), balls.end(), [&](const Ball &o) {
                return (o.c - b.c).norm() < o.R + b.R + 30.0;
            });

            if (clear) {
                balls.push_back(b);
                colour++;
            }
        }
    }

    void reset() {
        // the solver reads the balls, so any solve still running needs to finish first
        if (pending.valid()) {
            pending.get();
        }

        place_balls();
        run++;

        const Vector2d p0 = course.pts[0];
        body->state = trim.state(Vector3d(p0.x(), p0.y(), -h_c)); // trim faces north
        body->controls = body->trim;

        lat_ap.reset();
        lat_ap.clock.next = t;
        lon.reset();
        lon.clock.next = t;

        s_progress = 0.0;
        closest_pass = INF;
        trail.clear();

        have_plan = false;
        replan(true);
    }

    // the plane's (n, e, chi) with chi unwrapped to be near the plan's
    VectorXd measured_xi(const ControlState &cs, double chi_near) const {
        const double chi = cs.chi + 2.0 * std::numbers::pi *
                                        std::round((chi_near - cs.chi) / (2.0 * std::numbers::pi));
        return Vector3d(cs.pn, cs.pe, chi);
    }

    void replan(bool cold) {
        const ControlState cs = true_control_state(body->state, body->air);

        if (cold || !have_plan) {
            // straight ahead along the current course, no turning
            const VectorXd xi = measured_xi(cs, cs.chi);

            Z = {};
            for (int k = 0; k < s.N; k++) {
                VectorXd x = VectorXd::Zero(5);
                x.head<2>() =
                    xi.head<2>() + Va * k * s.dt() * Vector2d(std::cos(xi(2)), std::sin(xi(2)));
                x(2) = xi(2);

                Z.x.push_back(x);
                Z.u.push_back(VectorXd::Zero(1));
            }
            s.j_max = 30;
        } else {
            Z = SCvx::warm_start(ocp, s, t_plan, Z, 1);
            s.j_max = 4; // real time iteration-ish, the warm start is already close
        }

        const VectorXd xi_c = measured_xi(cs, Z.x[0](2));

        // V t_h further along the course from where the plane is now
        const double s_carrot = s_progress + Va * s.t_h;
        const Vector2d carrot = course.at(s_carrot);
        const double chi_carrot =
            course.heading(s_carrot) +
            2.0 * std::numbers::pi *
                std::round((xi_c(2) - course.heading(s_carrot)) / (2.0 * std::numbers::pi));

        s.xi_h = Vector3d(carrot.x(), carrot.y(), chi_carrot);

        // everything the solve needs is copied in, apart from the ocp (which only changes in reset)
        auto job = [this, settings = s, tc = t, xi_c, guess = Z]() {
            const auto t0 = std::chrono::steady_clock::now();
            Solved out{SCvx::solve(ocp, settings, tc, xi_c, guess), 0.0};
            out.ms =
                std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0)
                    .count();
            return out;
        };

        pending_t = t;
        if (async && !cold) {
            pending = std::async(std::launch::async, job);
        } else {
            adopt(job());
        }
    }

    // swaps the finished plan in
    void adopt(const Solved &done) {
        last = done.res;
        solve_ms = done.ms;
        Z = last.Z;
        t_plan = pending_t;
        have_plan = true;
    }

    // physics step
    void step(double dt) {
        if (pending.valid() &&
            pending.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
            adopt(pending.get());
        }

        if (have_plan && !pending.valid() && t >= t_plan + s.dt() - 1e-9) {
            replan(false);
        }

        const ControlState cs = true_control_state(body->state, body->air);

        // the lateral autopilot's loops, run here so the plan's turn rate can go straight into
        // the roll command as a feedforward, phi = atan(V r / g) (a coordinated turn)
        //
        if (lat_ap.clock.ready(t)) {
            const double Ts = lat_ap.clock.Ts;
            const double chi_ref = plan_course(Z, t_plan, s.dt(), t);
            const Vector2d p_ref = plan_position(Z, t_plan, s.dt(), ocp.V, t);
            const Vector2d right(-std::sin(chi_ref), std::cos(chi_ref));
            const double e_y = (Vector2d(cs.pn, cs.pe) - p_ref).dot(right);
            const double chi_c =
                wrap_angle(plan_course(Z, t_plan, s.dt(), t + lead) -
                           chi_inf * 2.0 / std::numbers::pi * std::atan(k_ct * e_y));
            const double phi_ff = std::atan(Va * plan_rate(Z, t_plan, s.dt(), t + lead) / 9.81);

            lat_ap.phi_c = lat_ap.course.pid.update(chi_c, cs.chi, Ts, phi_ff);
            body->controls.da = lat_ap.roll.update(lat_ap.phi_c, cs.phi, cs.p, Ts);

            if (lat_ap.has_rudder) {
                body->controls.dr = lat_ap.sideslip.update(cs.beta, Ts);
            }
        }
        lon.step(t, h_c, Va, cs.h, cs.Va, cs.theta, cs.q, body->controls);

        world.step(dt);
        t += dt;

        // progress only moves forward
        const Vector3d p = body->state.pos();
        double s_new;
        course.closest(p.head<2>(), &s_new);
        s_progress = std::max(s_progress, s_new);

        for (const Ball &b : balls) {
            const double d = (p - Vector3d(b.c.x(), b.c.y(), -h_c)).norm() - b.R;
            closest_pass = std::min(closest_pass, d);
            if (d < 0.0 && d > -dt * Va) {
                hits++;
            }
        }

        if (t >= next_trail) {
            trail.push_back(p);
            next_trail = t + 0.2;
        }

        if (s_progress >= course.length()) {
            reset();
        }
    }
};

// --- drawing ---
Vector3 rl(const Vector2d &ne, double h) {
    return rendering::ned_to_rl(Vector3d(ne.x(), ne.y(), -h));
}

const Color BALL_COLOURS[] = {RED, BLUE, GREEN, PURPLE, ORANGE, SKYBLUE, MAGENTA, GOLD, LIME, PINK};

void draw_world(const Mission &m) {
    for (size_t i = 0; i + 1 < m.course.pts.size(); i++) {
        DrawLine3D(rl(m.course.pts[i], 0.05), rl(m.course.pts[i + 1], 0.05), DARKGRAY);
        DrawLine3D(rl(m.course.pts[i], m.h_c), rl(m.course.pts[i + 1], m.h_c), LIGHTGRAY);
    }
    DrawCylinder(rl(m.course.pts.back(), 0.0), 0.5f, 0.5f, static_cast<float>(m.h_c), 8, DARKGREEN);

    // flown path
    for (size_t i = 0; i + 1 < m.trail.size(); i++) {
        DrawLine3D(rendering::ned_to_rl(m.trail[i]), rendering::ned_to_rl(m.trail[i + 1]),
                   DARKBLUE);
    }

    // the horizon, a straight line between each node
    // starts at the plane, since the first node or so is where it was when it last replanned
    // (up to a second behind it)
    const Vector3d pos = m.body->state.pos();
    const double h = -pos.z();

    Vector3 prev = rendering::ned_to_rl(pos);
    for (size_t k = 0; k < m.Z.x.size(); k++) {

        if (m.t_plan + static_cast<double>(k) * m.s.dt() <= m.t) {
            continue; // already past
        }

        const Vector3 node = rl(m.Z.x[k].head<2>(), h);
        DrawCylinderEx(prev, node, 0.08f, 0.08f, 6, ORANGE);
        DrawSphere(node, 0.3f, ORANGE);
        prev = node;
    }

    DrawSphereWires(rl(m.s.xi_h.head<2>(), h), 1.0f, 6, 6, DARKGREEN);
}

// semi transparent, so drawn last, back to front, without writing depth
void draw_balls(const Mission &m) {
    std::vector<const Ball *> order;
    for (const Ball &b : m.balls) {
        order.push_back(&b);
    }
    const Vector3 cam = rendering::CAMERA.position;
    auto dist = [&](const Ball *b) { return Vector3Distance(cam, rl(b->c, m.h_c)); };
    std::sort(order.begin(), order.end(),
              [&](const Ball *a, const Ball *b) { return dist(a) > dist(b); });

    rlDisableDepthMask();
    for (const Ball *b : order) {
        DrawSphereEx(rl(b->c, m.h_c), static_cast<float>(b->R), 24, 24,
                     Fade(BALL_COLOURS[b->colour % 10], 0.35f));
    }
    rlEnableDepthMask();
}

} // namespace

void run_scvx_scene() {
    rendering::init_camera({-10.0, -6.0, -6.0}, {0, 0, 0}, {0, 0, -1}, 50.0);
    Model wing = LoadModelFromMesh(GenMeshCube(1.4f, 0.05f, 0.3f));
    Model fin = LoadModelFromMesh(GenMeshCube(0.02f, 0.15f, 0.2f));

    Mission m;
    m.async = true;
    bool overview = false;
    double cam_psi = 0.0; // smoothed heading for the chase camera
    bool paused = false;

    const double dt_sim = m.body->sensors.params.Ts;
    double t_behind = 0.0;

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_R)) {
            m.reset();
        }
        if (IsKeyPressed(KEY_C)) {
            overview = !overview;
        }
        if (IsKeyPressed(KEY_T)) {
            m.world.atmosphere.turbulence = !m.world.atmosphere.turbulence;
        }
        if (IsKeyPressed(KEY_SPACE)) {
            paused = !paused;
        }

        if (!paused) {
            t_behind += std::min(static_cast<double>(GetFrameTime()), 0.25);
            while (t_behind >= dt_sim) {
                m.step(dt_sim);
                t_behind -= dt_sim;
            }
        }

        const Vector3d pos = m.body->state.pos();
        const Matrix3d R_b2i = rotation(m.body->state.quat());
        wing.transform = rendering::to_rl_transform(R_b2i, pos);
        fin.transform = rendering::to_rl_transform(R_b2i, pos + R_b2i * Vector3d(-0.1, 0, -0.075));

        if (overview) {
            const Vector3d mid(-25.0, 25.0, 0.0);
            rendering::init_camera(mid + Vector3d(0, 0, -650.0), mid, {1, 0, 0}, 60.0);
        } else {
            const double psi = quat_to_euler(m.body->state.quat()).z();
            cam_psi += std::min(1.0, 3.0 * GetFrameTime()) * wrap_angle(psi - cam_psi);
            const Vector3d fwd(std::cos(cam_psi), std::sin(cam_psi), 0.0);
            rendering::init_camera(pos - 14.0 * fwd + Vector3d(0, 0, -4.0), pos + 20.0 * fwd,
                                   {0, 0, -1}, 55.0);
        }

        rendering::begin(RAYWHITE, !overview);
        DrawModel(wing, {0, 0, 0}, 1.0f, GRAY);
        DrawModel(fin, {0, 0, 0}, 1.0f, DARKGRAY);
        draw_world(m);
        draw_balls(m);
        rendering::end_3d();

        const SCvx::Iteration &it = m.last.history.back();
        int row = 40;
        auto line = [&row](const char *text, Color c = DARKGRAY) {
            DrawText(text, 10, row, 18, c);
            row += 22;
        };
        line(TextFormat("scvx nmpc   N %d, horizon %.0f s, replanning every %.1f s", m.s.N, m.s.t_h,
                        m.s.dt()));
        line(TextFormat("last solve  %zu iterations  %.1f ms  defect %.1e  J_y %.1e",
                        m.last.history.size(), m.solve_ms, m.last.defect, it.J_y));
        line(TextFormat("closest pass %.1f m   hits %d   run %d", m.closest_pass, m.hits, m.run),
             m.hits > 0 ? RED : DARKGRAY);
        line(TextFormat("progress %.0f / %.0f m   turbulence %s", m.s_progress, m.course.length(),
                        m.world.atmosphere.turbulence ? "on" : "off"));
        DrawText("[r] new balls   [c] camera   [t] turbulence   [space] pause", 10,
                 GetScreenHeight() - 30, 16, DARKGRAY);
        rendering::end();
    }
    UnloadModel(wing);
    UnloadModel(fin);
}

} // namespace BSIM::scenes
