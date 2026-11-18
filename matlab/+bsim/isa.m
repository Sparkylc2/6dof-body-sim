function [T, p, rho, a] = isa(h)
%ISA international standard atmosphere at geometric altitude h [m]
%   [T, p, rho, a] = bsim.isa(h)
%
%   temperature [K], pressure [Pa], density [kg/m^3], speed of sound [m/s], each the same size as h
%   (valid to 32 km, held at the 32 km values above that)
%   note matlab's atmosisa takes geopotential altitude, this takes geometric (~13 m apart at 9 km)
[T, p, rho, a] = bsim_mex('isa', h);
end
