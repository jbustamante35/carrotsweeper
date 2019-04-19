function [mline, dsts] = generateMidline(msk, tCrds)
%% generateMidline: runs Nathan's algorithm for generating midlines from masks
% Description of Nathan's algorithm. Blah blah distance transform and
% mask gradient blah blah.
%
% A bit more detail about how the algorithms work.
%
% Usage:
%   mline = generateMidline(msk, tCrds)
%
% Input:
%   msk: binary mask image
%   tCrds: x-/y-coordinates representing the tip of the object
%
% Output:
%   mline: x-/y-coordinates defining the midline from the given mask
%

%% Constants
% Distance transform filtering
DSK   = 31;
FSPEC = 7;

% Tracking point at gradient
MAX_STEP = 50000;    % Original 50000
RHO      = 20;       % Original 20
RAD      = pi / 2;   % Original pi / 2
PDENSITY = [20 200]; % Original [20 200]
WSIGMA   = 0.3;      % Original 0.3

%% Euclidean distancs transform to compute distances to the edge of mask
gI       = double(bwdist(msk));
gI       = imfilter(gI, fspecial('gaussian', [DSK DSK], FSPEC));
[g1, g2] = gradient(gI);

%
img = double(bwdist(msk));
g1  = -g1;
g2  = -g2;

%
t1 = ba_interp2(g1, tCrds(1), tCrds(2));
t2 = ba_interp2(g2, tCrds(1), tCrds(2));

% Identify direction to begin tracing
N     = [t2 t1];
N     = -N / norm(N);
T     = [N(2) -N(1)];
iDirc = [T ; N];

%% Trace midline from tip coordinate along gradient
mline = trackFromPointAtGradient_carrot(img, tCrds, iDirc, ...
    MAX_STEP, RHO, RAD, PDENSITY, WSIGMA)';

% Get distances from midline to contour (hopefully)
dsts = getDim(impixel(gI, mline(:,1), mline(:,2)), 1) * 2;

end