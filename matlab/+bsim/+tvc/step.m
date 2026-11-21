function [s1, out, ov] = step(T, s0, delta_cmd, load, thrust, dt)
%STEP move the hydraulic TVC on by dt (its stiff hydraulics are sub stepped inside)
%   [s1, out, ov] = bsim.tvc.step(T, s0, delta_cmd, load, thrust, dt)
%
%   T          struct from bsim.tvc.params, or packed (faster in loops)
%   s0         12x1 [delta (2) | delta_dot (2) | P_L (4) | u_v (4)], zeros(12,1) is at rest centred
%   delta_cmd  2x1 gimbal angles the controller wants [rad]
%   load       2x1 the vehicle's load on the nozzle [N m] (the rocket terms' nozzle_load)
%   thrust     [N], for the thrust offset moment
%   dt         [s]
%
%   s1   12x1 state at the end
%   out  struct: delta, delta_dot, delta_ddot (average over the step), l, l_cmd, F (pistons)
%   ov   the same as a flat 18x1 vector
[s1, out, ov] = bsim_mex('tvc_step', T, s0, delta_cmd, load, thrust, dt);
end
