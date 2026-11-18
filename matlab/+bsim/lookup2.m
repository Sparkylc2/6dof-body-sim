function zq = lookup2(x, y, Z, xq, yq)
%LOOKUP2 the c++ 2d table lookup (bilinear, held flat past the edges)
%   zq = bsim.lookup2(x, y, Z, xq, yq)
%
%   Z is numel(x) by numel(y), so x goes down the rows and y along the columns
%   (that's the transpose of what interp2 wants, interp2(x, y, Z.', xq, yq) is the same inside)
zq = bsim_mex('lookup2', x, y, Z, xq, yq);
end
