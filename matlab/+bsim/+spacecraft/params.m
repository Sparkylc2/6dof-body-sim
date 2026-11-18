function [P, x0] = params(name)
%PARAMS parameters and starting state of a built in spacecraft
%   [P, x0] = bsim.spacecraft.params('cubesat')
%
%   P   struct: m_empty, J_empty (3x3), r_tank, R_tank, wheels (A, J_w, tau_max, h_max),
%       env (gravity, j2, gravity_gradient), thrusters (struct array: pos, dir, T_max, Isp)
%   x0  17x1 [r (3) | v (3) | q (4) | w (3) | h_wheels (3) | m], r and v inertial
[P, x0] = bsim_mex('spacecraft_params', name);
end
