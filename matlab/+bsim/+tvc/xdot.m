function [sd, out, ov] = xdot(T, s, delta_cmd, load, thrust)
%XDOT the TVC's 12 state derivatives (for analysis, simulate with bsim.tvc.step)
%   [sd, out, ov] = bsim.tvc.xdot(T, s, delta_cmd, load, thrust)
%   out has delta_ddot as the instantaneous acceleration, see bsim.tvc.step for the rest
[sd, out, ov] = bsim_mex('tvc_xdot', T, s, delta_cmd, load, thrust);
end
