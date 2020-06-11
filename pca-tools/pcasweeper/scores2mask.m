function [msk, prf] = scores2mask(scrs, evecs, mns)
%% scores2mask: generate averaged binary mask from set of PC scores
%
%
% Usage:
%   msk = scores2mask(scrs, evecs, mns)
%
% Input:
%   scrs: set or subset of PC scores
%   evecs: eigenvectors of the original dataset
%   mns: means of the original dataset
%
% Output:
%   msk: averaged binary mask representing the width profile
%   prf: the averaged width profile
%
WID = 300;

avg = mean(scrs, 1); % Mean by Columns
prf = pcaProject(avg, evecs, mns, 'scr2sim');
msk = profile2mask(prf, WID);

end
