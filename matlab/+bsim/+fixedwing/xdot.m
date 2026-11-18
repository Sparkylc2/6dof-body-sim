function [xd, F, M] = xdot(P, x, u, air)
%XDOT state derivative of a fixed wing
%   [xd, F, M] = bsim.fixedwing.xdot(P, x, u, air)
%
%   P    parameters (struct from bsim.fixedwing.params, or packed)
%   x    13x1 [pn pe pd | u v w | e0 e1 e2 e3 | p q r]   (u v w is the ground velocity, body)
%   u    4x1  [de da dr dt], where the surfaces actually are (rad), throttle 0-1
%   air  8x1  [rho a | wind_i (3) | gust_b (3)], from bsim.air or bsim.isa_air
%
%   xd   13x1 state derivative
%   F, M 3x1  body frame forces and moments it came from (the accelerometer wants F)
[xd, F, M] = bsim_mex('fixedwing_xdot', P, x, u, air);
end
