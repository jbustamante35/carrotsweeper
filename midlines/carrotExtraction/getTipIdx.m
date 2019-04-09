function [tCrds, tIdx] = getTipIdx(crv)
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

%% Constants for smoothing algorithm
SMOOTH = 10;  % Original 30

%% Run with default constant values
try
    % Smooth contour to remove noise (hairs, loops, etc)
%     bak = crv;
%     SMOOTHSPAN   = 0.7;
%     SMOOTHMETHOD = 'sgolay';
%     crv          = segSmooth(crv, SMOOTHSPAN, SMOOTHMETHOD);
    
    % Get maximum of computed curvature around contour
    curvature = cwtK(crv, SMOOTH);
    [~, tIdx] = max(curvature.K);
    tCrds     = crv(tIdx, :);
    
catch e
    % Try again with original constant values
    fprintf(2, 'Default to original constants...\n%s\n', e.message);
    KSNIP   = 50; % Original 50
    SMOOTH1 = 15; % Original 15
    SMOOTH2 = 30; % Original 30
    
    % Get minimum of inverse of initial curvature?
    oInit     = cwtK(crv, SMOOTH1);
    [~, tIdx] = min(oInit.K);
    
    % Get curvature around contour with finer smoothing
    oFine        = cwtK(crv, SMOOTH2);
    [~, fineIdx] = min((oFine.K(tIdx - KSNIP : tIdx + KSNIP)));
    
    % Compare tips from initial and finer smoothing
    tIdx  = tIdx + (fineIdx - KSNIP - 1);
    tCrds = crv(tIdx, :);
    
end

end