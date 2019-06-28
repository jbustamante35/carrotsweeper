function [tTip, tIdx, mDsk, mSmt, mTip, mCrv, hTip, rCnt] = tipFinderPlus(msk, xi, yi, smoothrange, diskrange, ncrds)
%% tipFinderPlus: 
%
%
% Usage:
%   [tTip, tIdx, mDsk, mSmt, mMsk, mCnt, mCrv , mIdx, mTip, hTip, iTip, rCnt] = ...
%       tipFinderPlus(msk, xi, yi, smoothrange, diskrange, ncrds)
%
% Input:
%
%
% Output:
%
%
%

%% Set Constants
CRDSTHRESH = 150;
CWTFLT     = 'closed';
PIX        = 200;
ALPHA      = 30;
NSMOOTH    = 3;
PIX2       = 40;
ALPHA2     = 1;

%% Get optimal parameters by iterating through smooth and disk ranges
[mDsk, mSmt] = ...
    iterativeSmoothing(msk, xi, yi, smoothrange, diskrange, ncrds, CWTFLT);

[mTip, mIdx, mCnt, mCrv] = ...
    getMaxParameters(msk, mDsk, mSmt, ncrds, CWTFLT, CRDSTHRESH);

%% Fix Initial guessed tip by smoothing around half-way region
% Get distances between coordinate indices then get pix / index scale
d   = diff(mCnt, 1, 1);
dL  = sum(d .* d, 2) .^ 0.5;
L   = cumsum([0 ; dL]);
SCL = L(end) / ncrds;

% Perform initial smoothing with large range and large smooth value
initCrv                  = mCrv;
initIdx                  = mIdx;
initCnt                  = [];
[tRgn, tCrv, tIdx, ~, tCnt] = ...
    tipSmoother(initCrv, initIdx, initCnt, mSmt, PIX, SCL, ALPHA, msk);

fprintf('|%d', tIdx);

% Decrease PIX range and ALPHA size to refine tip finder
for n = 1 : NSMOOTH
    [tRgn, tCrv, tIdx, ~, tCnt] = ...
        tipSmoother(tCrv, tIdx, tCnt, mSmt, PIX2, SCL, ALPHA2, []);
    
    fprintf('|%d|', tIdx);
end

%% Misc Data
hTip = tCnt(mIdx,:);        % Optimal tip index in original contour
rCnt = tCnt(tRgn > 0, :);   % Region to search for refined tip
tTip = tCnt(tIdx, :);       % Refined tip coordinates

end

function [mDsk, mSth] = iterativeSmoothing(msk, xi, yi, SMOOTH, DISKSIZE, NCRDS, CWTFLT)
%% iterativeSmoothing:
%
%
% Usage:
%
%
% Input:
%
%
% Output:
%
%

%%
S = zeros([numel(SMOOTH) , numel(DISKSIZE)]);

%% Iterate through ranges for SMOOTH value and DISK size
for e1 = 1 : numel(DISKSIZE)
    % Smooth binary mask and extract contour
    tmpmsk = imopen(msk, strel('disk', DISKSIZE(e1), 0));
    tmpcnt = extractContour(tmpmsk, NCRDS);
    tmpcnt.ReindexCoordinates('alt');
    tmpcnt = tmpcnt.NormalizedOutline;
    
    for e2 = 1 : numel(SMOOTH)
        % Compute curvature probabilities
        [~, crv] = cwtK(tmpcnt, SMOOTH(e2), CWTFLT);        
        ci       = interp1(xi, yi, crv, 'linear', 0);
        ci       = -log(ci);
        
        % Replace Inf values with 0 and sum probabilities
        ci(isinf(ci)) = 0;
        S(e1,e2)      = sum(ci);
    end
end

%% Get results using optimized parameters
% Get max of smatrix for carrot
[~ , maxS]    = max(S(:));
[mRow , mCol] = ind2sub(size(S), maxS);
mDsk          = DISKSIZE(mRow);
mSth          = SMOOTH(mCol);

end

function [mTip, mIdx, mCnt, mCrv] = getMaxParameters(msk, mDsk, mSmt, ncrds, cwtflt, crdsthresh)
%% getMaxParameters:
%
%
% Usage:
%
%
% Input:
%
%
% Output:
%
%

%% Perform operations using optimal parameters
% Open mask and extract contour
mMsk = imopen(msk, strel('disk', mDsk, 0));
mCnt = extractContour(mMsk, ncrds);
mCnt.ReindexCoordinates('alt');
mCnt = mCnt.NormalizedOutline;

% Smooth contoure and compute curvature
[~, mCrv] = cwtK(mCnt, mSmt, cwtflt);

% Filter out x-coordinates > threshold
crsdflt = mCnt(:,1) >= crdsthresh;
mCrv    = mCrv .* crsdflt;

% Get index and coordinates from point of maximum curvature
[~ , mIdx] = max(mCrv);
mTip       = mCnt(mIdx, :);
end
