function air = air(rho, a, wind_i, gust_b)
%AIR the 8x1 air vector [rho a | wind_i | gust_b] from given values
%   air = bsim.air()                         the c++ defaults (the book's rho, sea level a)
%   air = bsim.air(rho, a, wind_i, gust_b)   anything left out uses the defaults / zero wind
%
%   use bsim.isa_air to get rho and a from altitude instead
if nargin < 1 || isempty(rho), rho = 1.2682; end
if nargin < 2 || isempty(a), a = 340.29399; end
if nargin < 3, wind_i = zeros(3, 1); end
if nargin < 4, gust_b = zeros(3, 1); end
air = [rho; a; wind_i(:); gust_b(:)];
end
