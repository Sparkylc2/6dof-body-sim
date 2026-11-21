function T = unpack(p)
%UNPACK a packed TVC back into the struct (see bsim.tvc.params)
%   T = bsim.tvc.unpack(p)
T = bsim_mex('tvc_unpack', p);
end
