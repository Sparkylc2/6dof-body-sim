function air = isa_air(h, wind_i, gust_b)
%ISA_AIR the 8x1 air vector at altitude h [m] using the standard atmosphere
%   air = bsim.isa_air(h)                    no wind
%   air = bsim.isa_air(h, wind_i, gust_b)    steady wind (inertial) and gust (body), both 3x1
if nargin < 2, wind_i = zeros(3, 1); end
if nargin < 3, gust_b = zeros(3, 1); end
air = bsim_mex('isa_air', h, wind_i, gust_b);
end
