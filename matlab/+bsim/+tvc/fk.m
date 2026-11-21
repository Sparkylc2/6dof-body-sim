function [delta, residual] = fk(T, l, delta0)
%FK gimbal angles from four measured piston lengths (forward kinematics, gauss newton)
%   [delta, residual] = bsim.tvc.fk(T, l, delta0)
%   residual is the leftover length misfit [m]: big means a sensor disagrees with the others
if nargin < 3, delta0 = [0; 0]; end
[delta, residual] = bsim_mex('tvc_fk', T, l, delta0);
end
