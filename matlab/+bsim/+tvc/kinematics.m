function [l, J] = kinematics(T, delta)
%KINEMATICS piston lengths and the jacobian at gimbal angles delta (inverse kinematics)
%   [l, J] = bsim.tvc.kinematics(T, delta)
%   l  4x1 piston lengths [m],  J  4x2 with l_dot = J delta_dot (and gimbal torque = J' F)
[l, J] = bsim_mex('tvc_kinematics', T, delta);
end
