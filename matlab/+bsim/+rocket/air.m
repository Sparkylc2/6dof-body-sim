function a = air(h, wind_i, gust_b)
%AIR the rocket's 9x1 air vector at geometric altitude h [m] from the standard atmosphere
%   a = bsim.rocket.air(h)                    no wind
%   a = bsim.rocket.air(h, wind_i, gust_b)    steady wind (NED) and gust (body), both 3x1
%
%   a = [rho a p | wind_i | gust_b]
if nargin < 2, wind_i = zeros(3, 1); end
if nargin < 3, gust_b = zeros(3, 1); end
a = bsim_mex('isa_rocket_air', h, wind_i, gust_b);
end
