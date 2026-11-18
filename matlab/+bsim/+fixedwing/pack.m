function p = pack(P)
%PACK the parameter struct as the 47 long vector the c api takes
%   p = bsim.fixedwing.pack(P)
%
%   the packed vector skips the struct conversion on every call, so use it inside loops (scvx,
%   ode45) and as the parameter in simulink. errors on missing or misspelt fields
p = bsim_mex('fixedwing_pack', P);
end
