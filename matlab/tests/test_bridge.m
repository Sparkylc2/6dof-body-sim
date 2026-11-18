function tests = test_bridge
%TEST_BRIDGE checks the c++ bridge from the matlab side
%   run with: runtests('matlab/tests')   (after setup)
tests = functiontests(localfunctions);
end

function setupOnce(~)
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setup.m'));
end

% a zagi in roughly level flight at 20 m/s
function [P, x, u, air] = level_zagi()
P = bsim.fixedwing.params('zagi');
x = [0; 0; -100; 20; 0; 1; 1; 0; 0; 0; 0; 0; 0];
u = [-0.1; 0; 0; 0.6];
air = bsim.air();
end

% --- fixed wing ---
function test_pack_round_trip(tc)
P = bsim.fixedwing.params('zagi');
p = bsim.fixedwing.pack(P);
verifySize(tc, p, [47 1]);
verifyEqual(tc, bsim.fixedwing.unpack(p), P);
end

function test_unknown_field_is_an_error(tc)
P = bsim.fixedwing.params('zagi');
P.CL_alhpa = 1;
verifyError(tc, @() bsim.fixedwing.pack(P), 'bsim:error');
end

function test_struct_and_packed_agree(tc)
[P, x, u, air] = level_zagi();
verifyEqual(tc, bsim.fixedwing.xdot(P, x, u, air), ...
    bsim.fixedwing.xdot(bsim.fixedwing.pack(P), x, u, air));
end

% xdot is the wrench pushed through the equations of motion
function test_xdot_matches_wrench_and_eom(tc)
[P, x, u, air] = level_zagi();
[xd, F, M] = bsim.fixedwing.xdot(P, x, u, air);
[F2, M2] = bsim.fixedwing.wrench(P, x, u, air);
verifyEqual(tc, [F; M], [F2; M2]);
verifyEqual(tc, xd, bsim.aircraft_xdot(P.m, P.J, x, F, M), 'AbsTol', 1e-12);
% wings level heading north, so the position moves at [u; v; w] unrotated
verifyEqual(tc, xd(1:3), x(4:6), 'AbsTol', 1e-12);
end

% gravity only (no air), so the only acceleration is g down
function test_free_fall(tc)
P = bsim.fixedwing.params('zagi');
x = [0; 0; -100; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0];
xd = bsim.fixedwing.xdot(P, x, [0; 0; 0; 0], bsim.air(0));
verifyEqual(tc, xd(4:6), [0; 0; 9.81], 'AbsTol', 1e-12);
end

function test_jac_matches_matlab_finite_differences(tc)
[P, x, u, air] = level_zagi();
[A, B] = bsim.fixedwing.jac(P, x, u, air);
f = @(x, u) bsim.fixedwing.xdot(P, x, u, air);
h = 1e-5;
A_fd = zeros(13);
for i = 1:13
    e = zeros(13, 1);
    e(i) = h;
    A_fd(:, i) = (f(x + e, u) - f(x - e, u)) / (2 * h);
end
B_fd = zeros(13, 4);
for i = 1:4
    e = zeros(4, 1);
    e(i) = h;
    B_fd(:, i) = (f(x, u + e) - f(x, u - e)) / (2 * h);
end
verifyEqual(tc, A, A_fd, 'AbsTol', 1e-5);
verifyEqual(tc, B, B_fd, 'AbsTol', 1e-5);
end

% --- spacecraft ---
function test_spacecraft_round_trip_and_xdot(tc)
[P, x0] = bsim.spacecraft.params('cubesat');
verifySize(tc, x0, [17 1]);
p = bsim.spacecraft.pack(P);
verifySize(tc, p, [30 + 8 * numel(P.thrusters), 1]);
verifyEqual(tc, bsim.spacecraft.unpack(p), P);

% coasting with the wheels off, so only gravity: the position derivative is the velocity and the
% acceleration points at the earth
xd = bsim.spacecraft.xdot(P, x0, zeros(numel(P.thrusters), 1), zeros(3, 1));
verifyEqual(tc, xd(1:3), x0(4:6));
verifyLessThan(tc, dot(xd(4:6), x0(1:3)), 0);
verifyEqual(tc, xd(17), 0);
end

% --- atmosphere ---
% against matlab's own standard atmosphere (aerospace toolbox), which takes geopotential altitude
% --- rocket ---
function test_rocket_round_trip(tc)
P = bsim.rocket.params('example');
p = bsim.rocket.pack(P);
verifyEqual(tc, bsim.rocket.unpack(p), P);
verifyEqual(tc, bsim.rocket.pack(p), p); % packed goes straight through
end

function test_rocket_bad_input_is_an_error(tc)
[P, x0] = bsim.rocket.params('example');
u = [0; 0; 1; Inf; 0; zeros(P.fins.n, 1)];
air = bsim.rocket.air(0);

Q = P;
Q.engine.T_vacmax = 1;
verifyError(tc, @() bsim.rocket.pack(Q), 'bsim:error');
Q = P;
Q.drag.friction(2, 1) = Q.drag.friction(1, 1); % breakpoints not increasing
verifyError(tc, @() bsim.rocket.pack(Q), 'bsim:error');

p = bsim.rocket.pack(P);
verifyError(tc, @() bsim.rocket.xdot(p(1:end-1), x0, u, air, 0), 'bsim:error');
verifyError(tc, @() bsim.rocket.xdot(p, x0, u(1:end-1), air, 0), 'bsim:error');
end

% full throttle on the pad: not moving yet so there's no aero, just thrust straight up and gravity
function test_rocket_on_the_pad(tc)
[P, x0] = bsim.rocket.params('example');
u = [0; 0; 1; Inf; 0; zeros(P.fins.n, 1)];
air = bsim.rocket.air(0);
[xd, terms, tv] = bsim.rocket.xdot(P, x0, u, air, 0);

verifyEqual(tc, bsim.rocket.xdot(bsim.rocket.pack(P), x0, u, air, 0), xd);
verifyEqual(tc, terms.m, x0(14));
verifyEqual(tc, terms.V, 0);
verifyEqual(tc, terms.F_aero, zeros(3, 1));

% vacuum thrust less the exit pressure term, and the mass flow from the vacuum isp
T = P.engine.T_vac_max - air(3) * P.engine.A_e;
verifyEqual(tc, terms.T, T, 'RelTol', 1e-12);
verifyEqual(tc, xd(14), -P.engine.T_vac_max / (P.engine.Isp_vac * 9.80665), 'RelTol', 1e-9);

% nose up, so the thrust is along -down
verifyEqual(tc, terms.F_body, [T; 0; 0], 'AbsTol', 1e-9);
verifyEqual(tc, xd(4:6), terms.g_i - [0; 0; T / x0(14)], 'AbsTol', 1e-9);
verifyEqual(tc, xd(1:3), zeros(3, 1));
verifyLessThan(tc, xd(6), 0); % going up

verifyEqual(tc, tv, [terms.m; terms.x_cm; terms.Ix; terms.It; terms.T; terms.m_dot_e; ...
    terms.V; terms.mach; terms.qbar; terms.alpha_t; terms.F_body; terms.M_body; ...
    terms.F_aero; terms.M_aero; terms.g_i; terms.n_T; terms.rho_e; terms.v_r]);
verifyEqual(tc, terms.n_T, [1; 0; 0]);
end

% the viewer draws a frame from the packed P without complaining
function test_rocket_scene_draws(tc)
[P, x0] = bsim.rocket.params('example');
p = bsim.rocket.pack(P);
u = [deg2rad(3); 0; 1; Inf; 0; zeros(P.fins.n, 1)];
[~, terms] = bsim.rocket.xdot(p, x0, u, bsim.rocket.air(0), 0);
s = bsim.rocket.Scene(p, 'Visible', 'off', 'GimbalExaggeration', 2);
tc.addTeardown(@() close(s.Fig));
s.update(0, x0, terms);
s.update(0.1, x0, terms);
verifyTrue(tc, isvalid(s.Fig));
end

function test_rocket_air_layout(tc)
a = bsim.rocket.air(9000, [1; 2; 3], [4; 5; 6]);
[~, p, rho, s] = bsim.isa(9000);
verifyEqual(tc, a, [rho; s; p; 1; 2; 3; 4; 5; 6]);
end

function test_isa_matches_atmosisa(tc)
h = [0 1000 5000 9000 11500 15000 25000];
[T, p, rho, a] = bsim.isa(h);
Re = 6356766;
H = Re * h ./ (Re + h);
[T2, a2, p2, rho2] = atmosisa(H);
verifyEqual(tc, T, T2, 'RelTol', 1e-9);
verifyEqual(tc, a, a2, 'RelTol', 1e-7); % atmosisa's gas constant is a hair different
verifyEqual(tc, p, p2, 'RelTol', 1e-6);
verifyEqual(tc, rho, rho2, 'RelTol', 1e-6);
end

function test_isa_air_layout(tc)
air = bsim.isa_air(9000, [1; 2; 3], [4; 5; 6]);
[~, ~, rho, a] = bsim.isa(9000);
verifyEqual(tc, air, [rho; a; 1; 2; 3; 4; 5; 6]);
end

% --- tables ---
function test_lookup1_matches_interp1(tc)
x = [0 0.8 1.0 1.2 3.0];
y = [0.45 0.47 0.70 0.60 0.40];
xq = linspace(-1, 4, 51);
% interp1 with the ends clamped is the same thing
expected = interp1(x, y, min(max(xq, x(1)), x(end)));
verifyEqual(tc, bsim.lookup1(x, y, xq), expected, 'AbsTol', 1e-14);
end

function test_lookup2_matches_interp2(tc)
x = [0 0.5 0.9 1.1 2];
y = [0 0.05 0.1 0.2];
Z = x.' .^ 2 + 3 * y + x.' .* y;  % numel(x) by numel(y)
[XQ, YQ] = meshgrid(linspace(0, 2, 17), linspace(0, 0.2, 9));
verifyEqual(tc, bsim.lookup2(x, y, Z, XQ, YQ), interp2(x, y, Z.', XQ, YQ), 'AbsTol', 1e-13);
% past the edges it holds the edge
verifyEqual(tc, bsim.lookup2(x, y, Z, 5, -1), Z(end, 1), 'AbsTol', 1e-14);
end
