// src/scenes/sindy_scene.cpp
#include "SINDy/get_training_data.hpp"
#include "SINDy/model.hpp"
#include "body.hpp"
#include "frames.hpp"
#include "render/render.hpp"
#include "scenes/scenes.hpp"
#include "sensors.hpp" // DEG
#include "trim.hpp"
#include "world.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <raylib.h>
#include <vector>

namespace BSIM::scenes {
using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::MatrixXd;
using Eigen::Vector4d;
using Eigen::VectorXd;
using SINDy::SindyModel;

VectorXd rk4_step(const SINDy::SindyModel &model, const VectorXd &x, double dt) {
    const VectorXd k1 = model.f(x);
    const VectorXd k2 = model.f(x + dt / 2.0 * k1);
    const VectorXd k3 = model.f(x + dt / 2.0 * k2);
    const VectorXd k4 = model.f(x + dt * k3);
    return x + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

namespace {

constexpr double DT = 0.01;      // sim and rollout step (same as the sensor Ts)
constexpr double T_TRAIN = 30.0; // length of the sim
constexpr int N_TRAIN = 15;      // number of trajectories
constexpr double LAMBDA = 0.01;  // stlsq threshold (js the same as the tut)
constexpr double VA_TRIM = 18.0;
constexpr double H0 = 100.0; // initial height

// --- training ---
// one model is fit on the exact derivatives,
// and one on the central differences
struct Models {
    SindyModel exact, fd;
};

Models fit_models(const LongitudinalTrim &trim) {
    std::mt19937 rng(0);
    std::uniform_real_distribution<double> U(-1.0, 1.0);

    std::vector<SINDy::Trajectory<4>> trajs;
    int rows = 0;
    for (int i = 0; i < N_TRAIN; i++) {
        // random kick in all four states
        const double du = 3.0 * U(rng);        // m/s
        const double dw = 1.0 * U(rng);        // m/s
        const double q = 0.5 * U(rng);         // rad/s
        const double dth = 5.0 * DEG * U(rng); // degrees

        trajs.push_back(SINDy::simulate_longitudinal(trim, Vector4d(du, dw, q, dth), T_TRAIN, DT));
        rows += trajs.back().X.rows();
    }

    // all the different trajectories get stacked into these
    MatrixXd X(rows, 4);
    MatrixXd Xdot(rows, 4);
    MatrixXd Xdot_true(rows, 4);

    int r = 0;
    for (const auto &t : trajs) {
        const int m = static_cast<int>(t.X.rows());

        X.middleRows(r, m) = t.X;
        Xdot.middleRows(r, m) = t.Xdot;
        Xdot_true.middleRows(r, m) = t.Xdot_true;

        r += m;
    }

    Models out;

    out.exact.lib.state_names = {"du", "dw", "q", "dth"};
    out.fd.lib.state_names = out.exact.lib.state_names;

    out.exact.fit(X, Xdot_true, LAMBDA);
    out.fd.fit(X, Xdot, LAMBDA);
    return out;
}

// --- checks ---
// finite difference jacobian of the true longitudinal dynamics at trim
Matrix4d true_jacobian(const LongitudinalTrim &trim) {

    FixedWing b = generate_zagi_flying_wing();
    b.controls = trim.controls;

    const Vector4d x0 = SINDy::to_longitudinal(trim.state());

    auto f = [&b](const Vector4d &x) {
        b.state = SINDy::from_longitudinal(x);
        return SINDy::longitudinal_derivative(b);
    };

    const double h = 1e-6;
    const Vector4d f0 = f(x0);

    Matrix4d A;
    for (int j = 0; j < 4; j++) {
        Vector4d xp = x0;
        xp(j) += h;
        A.col(j) = (f(xp) - f0) / h;
    }
    return A;
}

// the model's jacobian at trim will just be the linear terms (approximately)
Matrix4d model_jacobian(const SindyModel &model) { return model.Xi.middleRows(1, 4).transpose(); }

// prints each mode (complex pairs only with +imag) with the natural frequency and damping
void print_modes(const char *label, const Matrix4d &A) {
    const Eigen::EigenSolver<Matrix4d> es(A);

    std::printf("  %-14s", label);

    for (int i = 0; i < 4; i++) {
        const std::complex<double> l = es.eigenvalues()(i);

        if (l.imag() < 0.0) {
            continue;
        }

        const double wn = std::abs(l);
        std::printf("[%+.3f %+.3fi    wn %.3f     zeta %.3f]", l.real(), l.imag(), wn,
                    -l.real() / wn);
    }
    std::printf("\n");
}

// prints rms error per state and the period from the same perturbation at the end of 40 s
void rms_period(const char *label, const SindyModel &model, const LongitudinalTrim &trim,
                const Vector4d &dx0) {

    const SINDy::Trajectory<4> traj = SINDy::simulate_longitudinal(trim, dx0, 40.0, DT);
    const int m = static_cast<int>(traj.X.rows());

    // row 0 of the sim is at t = dt
    MatrixXd Xs(m, 4);
    VectorXd x = traj.X.row(0).transpose();

    // steps the model
    for (int k = 0; k < m; k++) {
        Xs.row(k) = x.transpose();
        x = rk4_step(model, x, DT);
    }

    const Eigen::RowVectorXd rms =
        (Xs - traj.X).colwise().norm() / std::sqrt(static_cast<double>(m));

    std::printf("%-14s rms u %.4f m/s     w %.4f m/s  q %.4f rad/s    th %.4f deg", label, rms(0),
                rms(1), rms(2), rms(3) / DEG);
}

} // namespace

void run_sindy_scene() {
    BeginDrawing();
    ClearBackground(RAYWHITE);
    DrawText("generating and fitting", 20, 20, 20, DARKGRAY);
    EndDrawing();

    // --- trim and train ---
    const LongitudinalTrim trim = trim_longitudinal(generate_zagi_flying_wing(), VA_TRIM);

    if (!trim.converged) {
        TraceLog(LOG_WARNING, "trim didn't converge, residual %g", trim.residual);
    }

    const Vector4d x0 = SINDy::to_longitudinal(trim.state());
    const Models models = fit_models(trim);

    // --- startup report ---
    std::printf("\n-- SINDy, zagi at %.0f m/s, %d trajectories of %.0f s, lambda %g --- \n",
                VA_TRIM, N_TRAIN, T_TRAIN, LAMBDA);
    std::printf("\n     fit on exact derivatives:\n");
    models.exact.print();
    std::printf("\n     fit on central differences\n");
    models.fd.print();

    std::printf("\n-- modes at trim (short period, phugoid)\n");
    print_modes("true sim:", true_jacobian(trim));
    print_modes("sindy exact:", model_jacobian(models.exact));
    print_modes("sindy fd:", model_jacobian(models.fd));

    std::printf("\n-- result after 40 s (from a [+2 m/s, 0, 0, 0] perturb)\n");
    rms_period("sindy exact:", models.exact, trim, Vector4d(2.0, 0, 0, 0));
    std::printf("\n");
    rms_period("sindy fd:", models.fd, trim, Vector4d(2.0, 0, 0, 0));
    std::printf("\n\n");

    // --- drawing ---
    Mesh mesh = GenMeshCube(1.4f, 0.05f, 0.3f);
    Mesh mesh_fin = GenMeshCube(0.02f, 0.15f, 0.2f);
    Model model = LoadModelFromMesh(mesh);
    Model model_fin = LoadModelFromMesh(mesh_fin);

    auto draw_plane = [&](const Matrix3d &R, const Vector3d &pos, Color wing, Color fin) {
        model.transform = rendering::to_rl_transform(R, pos);
        model_fin.transform = rendering::to_rl_transform(R, pos + R * Vector3d(-0.1, 0, -0.075));
        DrawModel(model, {0, 0, 0}, 1.0f, wing);
        DrawModel(model_fin, {0, 0, 0}, 1.0f, fin);
    };

    World world;
    world.atmosphere.turbulence = false;
    auto plane = std::make_unique<FixedWing>(generate_zagi_flying_wing());
    FixedWing &body = *plane;
    world.bodies.bodies.push_back(std::move(plane));
    body.trim = trim.controls;

    // --- plane copy ---
    // exact derivative model +3 in y , central difference model +6 in y
    struct Ghost {
        const SindyModel *model;
        Vector3d offset;
        Color wing, fin;
        Vector4d xs = Vector4d::Zero(); // [du, dw, q, dth]
        Vector3d pos = Vector3d::Zero();
    };

    Ghost ghosts[] = {
        {&models.exact, Vector3d(0.0, 3.0, 0.0), Fade(BLUE, 0.7f), Fade(DARKBLUE, 0.7f)},
        {&models.fd, Vector3d(0.0, 6.0, 0.0), Fade(DARKGREEN, 0.7f), Fade(DARKGREEN, 0.7f)},
    };

    // perturbations indexed by keys 1-4 so u can play around a lil
    struct Kick {
        const char *name;
        Vector4d dx;
    };
    const Kick kicks[] = {
        {"+2 m/s u", Vector4d(2.0, 0, 0, 0)},
        {"+3 deg theta", Vector4d(0, 0, 0, 3.0 * DEG)},
        {"+0.3 rad/s q", Vector4d(0, 0, 0.3, 0)},
        {"+8 m/s, +10 deg", Vector4d(8.0, 0, 0, 10.0 * DEG)},
    };

    int kick = 0;
    bool paused = false;
    double t_behind = 0.0;
    auto reset = [&]() {
        const Vector4d dx = kicks[kick].dx;
        body.state = SINDy::from_longitudinal(x0 + dx);
        body.controls = trim.controls;
        for (Ghost &g : ghosts) {
            g.xs = dx;
            g.pos = body.state.pos() + g.offset;
        }
        t_behind = 0.0;
    };
    reset();

    while (!WindowShouldClose()) {
        for (int i = 0; i < 4; i++) {
            if (IsKeyPressed(KEY_ONE + i)) {
                kick = i;
                reset();
            }
        }
        if (IsKeyPressed(KEY_R)) {
            reset();
        }
        if (IsKeyPressed(KEY_SPACE)) {
            paused = !paused;
        }

        // --- steps both at the same fixed dt ---
        if (!paused) {
            const double frame_dt = GetFrameTime();
            t_behind += std::min(frame_dt, 0.25);
        }
        while (t_behind >= DT) {
            world.step(DT);

            for (Ghost &g : ghosts) {
                g.xs = rk4_step(*g.model, g.xs, DT);

                // position from the kinematics (eqn 3.14 with phi = psi = 0)
                const Vector4d x = x0 + g.xs;
                const double u = x(0), w = x(1), th = x(3);
                g.pos.x() += DT * (u * std::cos(th) + w * std::sin(th));
                g.pos.z() += DT * (-u * std::sin(th) + w * std::cos(th));
            }
            t_behind -= DT;
        }

        const Vector3d target(body.state.pos().x(), 3.0, -H0);
        rendering::init_camera(target + Vector3d(-18.0, -22.0, -5.0), target, {0, 0, -1}, 45.0);

        rendering::begin();
        draw_plane(rotation(body.state.quat()), body.state.pos(), GRAY, DARKGRAY);
        for (const Ghost &g : ghosts) {
            draw_plane(rotation(euler_to_quat(0.0, x0(3) + g.xs(3), 0.0)), g.pos, g.wing, g.fin);
        }
        rendering::end_3d();

        DrawText("grey: sim", 10, 10, 14, DARKGRAY);
        DrawText("blue: sindy exact (3 m)", 90, 10, 14, BLUE);
        DrawText("green: sindy fd (6 m)", 290, 10, 14, DARKGREEN);
        rendering::end();
    }

    UnloadModel(model);
    UnloadModel(model_fin);
}

} // namespace BSIM::scenes
