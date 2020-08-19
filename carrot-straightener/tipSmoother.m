function [sVec, sCrv, sIdx, sTip, iCnt] = tipSmoother(iCrv, iTip, iCnt, optS, PIX, SCL, ALPHA, init)
%% tipSmoother: windowed smoothing operation to find tip location
%
%
% Usage:
%   [sVec, sCrv, sIdx, sTip, iCnt] = tipSmoother( ...
%               iCrv, iTip, iCnt, optS, PIX, SCL, ALPHA)
%
% Input:
%   iCrv: array of curvatures from closed contour
%   iTip: tip coordinates to draw a window from
%   iCnt: contour from original mask
%   optS: optimal smoothing value for computing curvature
%   PIX: number of pixel coordinates to generate window
%   SCL: scale factor to convert contour coordinates to image pixel distance
%   ALPHA: factor to scale smoothing parameter
%   init: original mask to initialize, if this is the first run
%
% Output:
%   sVec: window region to use for smoothing and locating the tip
%   sCrv: curvature values around window regions
%   sIdx: updated tip index after smoothing
%   sTip: updated tip coordiantes after smoothing
%   iCnt: return the original contour that was inputted

% Create initial curvature array and tip index
if ~isempty(init)
    % Extract contour from original mask
    iCnt  = extractContour(init, length(iCrv));
    iCnt.ReindexCoordinates('alt');
    iCnt  = iCnt.NormalizedOutline;
end

% Create window region to smooth and compute curvatures
WIN                 = PIX / SCL;
sVec                = zeros(size(iCrv));
tipUp               = ceil(iTip - WIN);
tipDn               = ceil(iTip + WIN);
sVec(tipUp : tipDn) = 1;

% Smooth window and compute curvature to find tip index and coordinates
tmpOpt    = ALPHA * optS;
[~, sCrv] = cwtK(iCnt, tmpOpt, 'closed');
sCrv      = sCrv .* sVec;

if sum(sCrv < 0) ~= sum(sVec == 1)
    % If no negative curvatures
    [~ , sIdx] = max(sCrv);    
else
    % Get maximum of negative curvatures if entire region is negative
    fprintf('*');
    [nCrv, ~] = max(sCrv(sVec == 1));
    sIdx      = find(sCrv == nCrv);
end

sTip = iCnt(sIdx,:);

end
