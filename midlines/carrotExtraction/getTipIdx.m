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
    
    %% 
%     set(0, 'CurrentFigure', figure(3));
%     cla;clf;
%     
%     n    = 30;
%     srt  = fliplr(sort(curvature.K));
%     topN = srt(1:n);
%     plot(topN);
%     xlim([-10 10]);
% %     ylim([-0.2 10]);
% 
%     
%     topcrvs = crv(ismember(curvature.K , topN), :);
%     plt(crv, 'k.', 10);
%     hold on;
%     plt(topcrvs(1, :), 'rx', 5);
%     plt(topcrvs(11, :), 'bx', 5);
%     plt(topcrvs(21, :), 'mx', 5);
%     axis image;
%     axis ij;
    
    %% Try secondary smoothing
%     % Get maximum of computed curvature around contour
% %     SMOOTH2    = round(SMOOTH / 5);
%     SMOOTH2    = 3;
%     cIdx = 1;
%     c = 1;
%     
%     for c = 1 : 10 : length(topcrvs)        
%         curvature2(cIdx) = cwtK(topcrvs(c:(c + 9), :), SMOOTH2);
%         [~, tIdx2] = max(curvature2(cIdx).K);
%         tCrds2     = topcrvs(tIdx2, :);
%         cIdx = cIdx + 1;
%     end
    
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