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
KSNIP   = 163; % Original 50
SMOOTH1 = 5;   % Original 15
SMOOTH2 = 10;  % Original 30

%% Run with default constant values
try
%     %
%     oInit     = cwtK(crv, SMOOTH1);
%     [~, initIdx] = min(oInit.K);
%     
%     %
%     oFine        = cwtK(crv, SMOOTH2);    
%     
%     try
%         [~, fineIdx] = min((oFine.K(tIdx - KSNIP : tIdx + KSNIP)));
%     catch e
%         [~, fineIdx] = max(oFine.K);
%     end
%     
%     %
%     tIdx  = tIdx + (fineIdx - KSNIP - 1);
%     tCrds = crv(tIdx, :);
    
    %
    oFine     = cwtK(crv, SMOOTH2);    
    [~, tIdx] = max(oFine.K);
    tCrds     = crv(tIdx, :);
    
catch e
    fprintf(2, 'Default to original constants...\n%s\n', e.message);
    KSNIP   = 50;
    SMOOTH1 = 15;
    SMOOTH2 = 30;
    
    % Try again with original constant values
    oInit     = cwtK(crv, SMOOTH1);
    [~, tIdx] = min(oInit.K);
    
    %
    oFine        = cwtK(crv, SMOOTH2);
    [~, fineIdx] = min((oFine.K(tIdx - KSNIP : tIdx + KSNIP)));
    
    %
    tIdx  = tIdx + (fineIdx - KSNIP - 1);
    tCrds = crv(tIdx, :);
    
end

end