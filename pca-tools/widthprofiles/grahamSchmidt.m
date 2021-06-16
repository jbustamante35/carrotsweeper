function H = grahamSchmidt(bw, bl, scrs)
%% grahamSchmidt: orthonormalization method
%
%
% Usage:
%   H = grahamSchmidt(bw, bl, scrs)
%
% Input:
%
%
% Output:
%
%

%%
nbw = bw / norm(bw);
nbl = bl / norm(bl);

% Decompose NBW into NBL
% nu = nbw' - (nbl' * nbw)' * nbl'; % is this more correct?
nu  = nbw - (nbl' * nbw) * nbl;
nnu = nu / norm(nu);

% H   = S - (S * nbw) * nbw' - (S * nnu) * nnu';
% H   = scrP - (scrP * nbl) * nbl' - (scrP * nnu) * nnu';
H1 = scrs - (scrs * nbl) * nbl';
H  = H1 - (H1 *  nnu) * nnu';

end