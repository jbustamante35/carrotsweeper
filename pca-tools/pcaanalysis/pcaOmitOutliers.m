function [Q , remIdx] = pcaOmitOutliers(P, outlier_pct, pcdim)
%% pcaOmitOutliers: run PCA on dataset using outlier-removed eigenvectors
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
%   Q: pca data run on full dataset with outliers removed
%   remIdx: indices of remaining data
%

if nargin < 2; outlier_pct = [5 , 95]; end
if nargin < 3; pcdim       = 1;        end

%% Initial PCA, Remove Outliers, Re-Run PCA
% Get PCScores and DataName
S   = P.PCAScores;
pca = regexpi(P.DataName, 'pcaResults_', 'end');
pcb = regexpi(P.DataName, '_\dPCs',      'start');
dnm = P.DataName(pca + 1 : pcb - 1);

% Remove outliers and redo PCA
[~ , wIdx]  = rmoutliers(S(:,pcdim), 'percentiles', outlier_pct);
remIdx      = find(~wIdx);
wx          = P.InputData(remIdx,:);
W           = pcaAnalysis(wx, P.NumberOfPCs, 0, dnm);

Q.DataName    = W.DataName;
Q.InputData   = P.InputData;
Q.NumberOfPCs = W.NumberOfPCs;
Q.EigVecs     = W.EigVecs;
Q.MeanVals    = W.MeanVals;
Q.PCAScores   = pcaProject(Q.InputData, W.EigVecs, W.MeanVals, 'sim2scr');
Q.SimData     = pcaProject(Q.PCAScores, W.EigVecs, W.MeanVals, 'scr2sim');
end
