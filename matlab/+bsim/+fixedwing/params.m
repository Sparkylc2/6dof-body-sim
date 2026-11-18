function P = params(name)
%PARAMS parameters of one of the built in fixed wings, as a struct
%   P = bsim.fixedwing.params('zagi')    the zagi flying wing (appendix E)
%   P = bsim.fixedwing.params('sport')   the sport drone
%
%   P has m, J ([Jx; Jy; Jz; Jxz]) and every aero/propulsion coefficient, change any of them and
%   pass P straight to the other bsim.fixedwing functions
%   (pack it with bsim.fixedwing.pack first if it's going into a loop or simulink)
P = bsim_mex('fixedwing_params', name);
end
