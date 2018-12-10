function [E, lambda, gamma, isConvex] = lineIntersection(A,B,C,D)
% [E, lambda, gamma, isConvex] = lineIntersection(A,B,C,D)
%
% Given a line segment AB and another line segment CD, compute the point E
% where the lines intersect.
%
% INPUTS:
%   A = [2,n] = [Ax;Ay] = point in 2D space
%   B = [2,n] = [Bx;By] = point in 2D space
%   C = [2,n] = [Cx;Cy] = point in 2D space
%   D = [2,n] = [Dx;Dy] = point in 2D space
%
% OUTPUTS:
%   E = [2, n] = intersection of lines AB and CD
%   lambda = [1,n]
%       E = lambda*A + (1-lambda)*B
%   gamma = [1,n]
%       E = gamma*C + (1-gamma)*D
%   isConvex = is intersection on both lines?
%       isConvex = (0 <= lambda <= 1)  && (0 <= gamma <= 1)
%
% DERIVATION:
%   E1 = lambda*A + (1-lambda)*B
%   E2 = gamma*C + (1-gamma)*D
%   E1 == E2  --> linear system in [lambda; gamma] --> solve 
%
% IMPLEMENTATION:
%   F = B-D;
%   M = [(B-A), (C-D)]
%   Z = M\F;
%   lambda = Z(1);
%   gamma = Z(2);
%
F1 = B(1,:)-D(1,:);
F2 = B(2,:)-D(2,:);
M11 = B(1,:)-A(1,:);
M21 = B(2,:)-A(2,:);
M12 = C(1,:)-D(1,:);
M22 = C(2,:)-D(2,:);
deter = M11.*M22 - M12.*M21;
lambda = -(F2.*M12-F1.*M22)./deter;
gamma = (F2.*M11-F1.*M21)./deter;
E = ([1;1]*lambda).*A + ([1;1]*(1-lambda)).*B;
isConvex = (0 <= lambda & lambda <= 1)  & (0 <= gamma & gamma <= 1) ;
% E_check = ([1;1]*gamma).*C + ([1;1]*(1-gamma)).*D;  % Should match E
end