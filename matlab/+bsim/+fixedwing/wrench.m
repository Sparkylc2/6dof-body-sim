function [F, M] = wrench(P, x, u, air)
%WRENCH body frame forces and moments on a fixed wing (gravity, aero, propulsion)
%   [F, M] = bsim.fixedwing.wrench(P, x, u, air)
%
%   same inputs as bsim.fixedwing.xdot
[F, M] = bsim_mex('fixedwing_wrench', P, x, u, air);
end
