function P = unpack(p)
%UNPACK the packed spacecraft vector back into the struct
%   P = bsim.spacecraft.unpack(p)
P = bsim_mex('spacecraft_unpack', p);
end
