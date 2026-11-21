function [xd, terms, tv] = xdot(P, x, u, air, t)
%XDOT state derivative of the rocket
%   [xd, terms, tv] = bsim.rocket.xdot(P, x, u, air, t)
%
%   P      struct from bsim.rocket.params, or packed (packed is faster, use it in loops)
%   x      14x1 [pn pe pd | vn ve vd | e0 e1 e2 e3 | p q r | m] (cm pos and vel, NED)
%   u      (9 + fins.n)x1 [gimbal_p gimbal_y throttle t_ignition fins_deployed | fin angles |
%                          gimbal_p_dot gimbal_y_dot gimbal_p_ddot gimbal_y_ddot]
%          where the actuators actually are (rad), throttle 0-1, t_ignition is the sim time a
%          thrust curve motor was lit (Inf for not yet, throttleable engines ignore it), and how
%          the gimbal is moving (from its servos, for the swinging engine's dynamics). the last
%          four can be left off (5 + fins.n long), which treats the gimbal as standing still
%   air    9x1 [rho a p | wind_i (3, NED) | gust_b (3, body)], see bsim.rocket.air
%   t      sim time [s]
%
%   xd     14x1
%   terms  struct of what went into it: m, x_cm, Ix, It, T, m_dot_e, V, mach, qbar, alpha_t,
%          F_body, M_body (every non gravitational force / moment about the cm, body frame,
%          F_body / m is what an accelerometer reads), F_aero, M_aero, g_i (NED),
%          n_T (thrust direction, follows the gimbal even when the engine's off),
%          rho_e (nozzle exit from the cm), v_r (velocity relative to the air, body),
%          M_gimbal (every moment from the swinging engine, already in M_body),
%          hinge (torque the gimbal actuator has to supply, [pitch; yaw]),
%          nozzle_load (the vehicle's load on the nozzle, for bsim.tvc.step's load)
%   tv     the same as a flat 41x1 vector, in the order above (what the c api gives simulink)
[xd, terms, tv] = bsim_mex('rocket_xdot', P, x, u, air, t);
end
