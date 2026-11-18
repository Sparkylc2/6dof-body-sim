function [A, B] = jac(P, x, u, air)
%JAC linearisation of the fixed wing about (x, u)
%   [A, B] = bsim.fixedwing.jac(P, x, u, air)
%
%   A = d(xdot)/dx (13x13), B = d(xdot)/du (13x4), from central finite differences in c++
%   same inputs as bsim.fixedwing.xdot
[A, B] = bsim_mex('fixedwing_jac', P, x, u, air);
end
