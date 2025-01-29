function [tCrds , tIdx , debug_output] = getTipIdx(msk, VIS, DEBUG)
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
%   DEBUG: returns additional data for debugging
%
% Output:
%   tCrds: x-/y-coordinate corresponding to tip
%   tIdx: index along input curve corresponding to tip
%

if nargin < 2
    VIS          = 0;
    DEBUG        = 0;
    debug_output = [];
end

%% Set Constants
SMOOTHRANGE = 7 : 1 : 13;
DISKRANGE   = 7 : 1 : 13;
NCRDS       = 800;

%% Pull down trained curvature distribution
% Load tip curvature distribution
cshome = fileparts(which('getTipIdx'));
% tnfin  = sprintf('training/%s_CurveDistribution_%s_%dObjects_%dCurves.mat', ...
%     '190627', 'Initialize', 100, 70834);
tnfin  = sprintf('%s%s%s', 'training', filesep, 'trainingdata.mat');
DOUT   = load(sprintf('%s/%s', cshome, tnfin), 'DOUT');
OUT    = DOUT.DOUT;
xi     = OUT.xi;
yi     = OUT.yi;

try
    %% Use TipFindrX method to find optimally-refined tip
    % WTF I need to clean this up further    
    tic;
    if VIS; fprintf('Running TipFindrX...'); end
    
    [tCrds , tIdx , mDsk , mSmt , mTip , mCrv , hTip , rCnt] = ...
        tipFinderPlus(msk, xi, yi, SMOOTHRANGE, DISKRANGE, NCRDS, VIS);
    
    if DEBUG
        fnms = {'max_disk_size' , 'max_smooth_size' , 'optimal_tip_coord' , ...
            'optimal_contour' , 'optimal_tip_idx' , 'optimal_region'};
        debug_output = cell2struct({mDsk, mSmt, mTip, mCrv, hTip, rCnt}', fnms);
    end
    
    if VIS; fprintf('...Identified tip in %.02f sec\n', toc); end
    
catch e
    fprintf(2, 'Error with tip finder\n%s\n', e.getReport);
    tCrds = [];
    tIdx  = [];
end
end
