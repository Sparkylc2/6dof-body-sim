function p = pack(P)
%PACK the spacecraft struct as the packed vector the c api takes
%   p = bsim.spacecraft.pack(P)    (30 + 8 per thruster long)
p = bsim_mex('spacecraft_pack', P);
end
