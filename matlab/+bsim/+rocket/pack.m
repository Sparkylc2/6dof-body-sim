function p = pack(P)
%PACK the rocket struct as the packed vector the c api takes (and simulink wants)
%   p = bsim.rocket.pack(P)
%
%   the length depends on how many slices, table points etc there are, see the layout in
%   matlab/bridge/bsim_c.h. passing an already packed vector checks it and hands it back
p = bsim_mex('rocket_pack', P);
end
