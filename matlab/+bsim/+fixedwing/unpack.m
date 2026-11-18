function P = unpack(p)
%UNPACK the packed 47 long vector back into the parameter struct
%   P = bsim.fixedwing.unpack(p)
P = bsim_mex('fixedwing_unpack', p);
end
