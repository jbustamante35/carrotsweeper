function [Q, remIdx] = pcaOmitOutliers(P, outlier_pct, pcdim)
%% pcaOmitOutliers:
%
%
% Usage:
%
%
% Input:
%   P: PcaJB object to remove outliers
%   outlier_pct: percentage of outliers to omit (default: [5 , 94])
%   pcdim: PC score to check for outliers (default: 1)
%
% Output:
%   Q: new PcaJB object with outliers removed
%   remIdx: indices of remaining data
%

if nargin < 2
    outlier_pct = [5 , 95];
    pcdim       = 1;
end

%% Initial PCA, Remove Outliers, Re-Run PCA
% Get PCScores and DataName
S   = P.PCAScores;
pca = regexpi(P.DataName, 'pcaResults_', 'end');
pcb = regexpi(P.DataName, '_\dPCs',      'start');
dnm = P.DataName(pca + 1 : pcb - 1);

% Remove outliers and redo PCA
[~ , wIdx] = rmoutliers(S(:,pcdim), 'percentiles', outlier_pct);
remIdx     = find(~wIdx);
W          = P.InputData(remIdx,:);
Q          = pcaAnalysis(W, P.NumberOfPCs, 0, dnm, 0);

end
