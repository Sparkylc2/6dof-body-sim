function yq = lookup1(x, y, xq)
%LOOKUP1 the c++ 1d table lookup (linear, held flat past the ends)
%   yq = bsim.lookup1(x, y, xq)
%
%   the same thing the c++ models do with a Table1D, so checking a table here checks what the
%   model sees. x strictly increasing, xq any size
%   (matlab's interp1 gives NaN past the ends unless told otherwise, this holds the end values)
yq = bsim_mex('lookup1', x, y, xq);
end
