function xd = aircraft_xdot(m, J, x, F, M)
%AIRCRAFT_XDOT rigid body equations of motion for a symmetric aircraft
%   xd = bsim.aircraft_xdot(m, J, x, F, M)
%
%   for when the forces come from somewhere other than the fixed wing model
%   m    mass
%   J    4x1 [Jx; Jy; Jz; Jxz]
%   x    13x1 state (see bsim.fixedwing.xdot)
%   F, M 3x1 body frame forces and moments
xd = bsim_mex('aircraft_xdot', m, J, x, F, M);
end
