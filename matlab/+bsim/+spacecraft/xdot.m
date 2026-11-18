function xd = xdot(P, x, throttle, wheel_torque)
%XDOT state derivative of a spacecraft
%   xd = bsim.spacecraft.xdot(P, x, throttle, wheel_torque)
%
%   P             struct from bsim.spacecraft.params, or packed
%   x             17x1 state (see bsim.spacecraft.params)
%   throttle      one 0-1 entry per thruster
%   wheel_torque  3x1 commanded motor torque per wheel (the limits are applied inside)
xd = bsim_mex('spacecraft_xdot', P, x, throttle, wheel_torque);
end
