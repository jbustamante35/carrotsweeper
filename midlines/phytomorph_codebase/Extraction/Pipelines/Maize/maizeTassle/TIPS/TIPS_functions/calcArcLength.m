function [ arcLength ] = calcArcLength( S )
%CALCARCLENGTH Calculates arc length of a piecewise polynomial
%   S is a pp spline

Sder = fnder(S);
nSder = size(Sder.coefs, 1);
xt = Sder.coefs(1:2:nSder, :);
yt = Sder.coefs(2:2:nSder, :);
ds = (xt.^2 + yt.^2).^(.5);
ds = mkpp(1:(size(ds, 1) + 1), ds);
arcLength = diff(fnval(fnint(ds),[1 ds.pieces]));

end

