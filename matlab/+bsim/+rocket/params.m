function [P, x0, act] = params(name)
%PARAMS parameters, starting state and actuator numbers of a built in rocket
%   [P, x0, act] = bsim.rocket.params('example')
%
%   'example' is the Ember stand in from generate_example_rocket (craft.cpp)
%
%   P   struct, the same fields as RocketParams (force_models/rocket_model.hpp) plus env:
%       m_dry, x_dry, Ix_dry, It_dry, d, A_base, A_motor, C_dc, CD_flat
%       propellant  struct array: kind ('liquid' | 'solid'), share, x, R, rho, L
%       engine      kind ('throttleable' | 'thrust_curve'), thrust_curve ([t T] columns),
%                   c_eff, p_test, T_vac_max, Isp_vac, throttle_min, A_e, x_G, L_n
%       components  struct array: CN_alpha, x_cp
%       slices      n x 3 [x d dx]
%       drag        friction, pressure, parasitic, each [mach cd] columns ([] is zero)
%       fins        n, x, r, phi0, S, CN_alpha, CD0
%       env         gravity, aero, jet_damping, J_dot, coriolis (logical)
%   x0  14x1 [pn pe pd | vn ve vd | e0 e1 e2 e3 | p q r | m], on the pad, nose up, full
%   act the grid fin servo numbers (fin: kind, tau, wn, zeta, rate_max, pos_max), for building
%       them in simulink. the gimbal is the hydraulic TVC, see bsim.tvc.params('ember')
[P, x0, act] = bsim_mex('rocket_params', name);
end
