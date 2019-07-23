function [tCrds, tIdx] = getTipIdx(msk)
%% getTipIdx: identify point of highest curvature on a contour
% This function computes the curvature around a given contour and identifies the
% index along the curve representing the highest curvature. Parameters can be
% set to control the smoothing and step sizes for optimizing the algorithm.
%
% Usage:
%   [tIdx , tCrds] = getTipIdx(crv)
%
% Input:
%   crv: x-/y-coordinates of the contour
%
% Output:
%   tCrds: x-/y-coordinate corresponding to tip
%   tIdx: index along input curve corresponding to tip
%

%% Set Constants
SMOOTHRANGE = 7 : 1 : 13;
DISKRANGE   = 7 : 1 : 13;
NCRDS       = 800;

%% Pull down trained curvature distribution
% Load tip curvature distribution
cshome = fileparts(which('getTipIdx'));
tnfin  = sprintf('training/%s_CurveDistribution_%s_%dObjects_%dCurves.mat', ...
    '190627', 'Initialize', 100, 70834);
DOUT   = load(sprintf('%s/%s', cshome, tnfin), 'DOUT');
OUT    = DOUT.DOUT;
xi     = OUT.xi;
yi     = OUT.yi;

%% Use TipFindrX method to find optimally-refined tip
% WTF I need to clean this up further
tic;
fprintf('Running TipFindrX...');
[tCrds, tIdx] = tipFinderPlus(msk, xi, yi, SMOOTHRANGE, DISKRANGE, NCRDS);
fprintf('...Identified tip in %.02f sec\n', toc);

end
