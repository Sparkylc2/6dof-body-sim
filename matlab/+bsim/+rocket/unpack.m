function S = unpack(p)
%UNPACK a packed rocket back into the struct (see bsim.rocket.params)
%   S = bsim.rocket.unpack(p)
S = bsim_mex('rocket_unpack', p);
end
