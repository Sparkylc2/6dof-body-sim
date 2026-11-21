function p = pack(T)
%PACK the TVC struct as the 23x1 vector the c api takes (and simulink wants)
%   p = bsim.tvc.pack(T)
p = bsim_mex('tvc_pack', T);
end
