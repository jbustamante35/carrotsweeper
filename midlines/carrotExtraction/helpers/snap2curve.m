function [snp, snpIdx] = snap2curve(crd, crds)
%% snap2curve: snap coordinates to closest point along curve
% This function takes [m x n] coordinate positions, finds the index along [p x n] coordinate matrix,
% and returns an [m x n] matrix, where coordinates are replaced by nearest coordinates in crds. 
%
% Usage:
%   [snp, idx] = snap2curve(crd, crds)
% 
% Input:
%   crd: [m x n] matrix of coordinates near curve 
%   crds: [p x n] matrix of coordinates on curve to search along
%
% Output:
%   snp: [m x n] matrix of coordinates corresponding to nearest point from pts on curve crd
%   snpIdx: index in crds that was snapped to
%

%% Find indices corresponding to nearest distance from coordinates in pts
idx    = dsearchn(crds, crd);
snp    = crds(idx, :);
snpIdx = find(sum(crds == snp, 2) == 2);
end

