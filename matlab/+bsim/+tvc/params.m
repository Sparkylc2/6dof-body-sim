function T = params(name)
%PARAMS parameters of a built in hydraulic TVC (force_models/tvc_model.hpp)
%   T = bsim.tvc.params('ember')
%
%   the placeholder Ember TVC from the TVC actuator notes, with the example rocket's engine as the
%   nozzle. fields (SI, measured from the gimbal pivot in body axes):
%     geometry    x_a r_a (airframe pins), x_b r_b (nozzle pins, straight), phi0, stroke
%     hydraulics  A_p P_s Q_max beta V_t C_tp tau_v K_p
%     nozzle      m_n l_n I_na I_nt
%     loads       offset_p offset_y (thrust line miss [m]), friction (coulomb [N m]), damping
%     h           the internal integration step [s]
T = bsim_mex('tvc_params', name);
end
