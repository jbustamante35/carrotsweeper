function [W, wids, lngs] = prepWidthProfiles(D, nrmL, nrmW, INTERP)
%% prepWidthProfiles: prepare width profiles for PCA
%
% NOTE:
%   This function flips the widths left-to-right!
%
% Usage:
%   [W, wids, lngs] = prepWidthProfiles(D, nrmL, nrmW, INTERP)
%
% Input:
%   D: cell array of width profiles
%   nrmL: boolean to set normalization method for lengths
%   nrmW: boolean to set normalization method for widths
%   I: size to interpolate length to [defaults to longest in dataset]
%
% Output:
%   W: width profiles stacked using requested method
%   wids: original widths for data regressor of Graham Schmidt analysis
%   lngs: original lengths for data regressor of Graham Schmidt analysis
%

%% Handle lengths
if nrmL
    % Set interpolation length to longest in dataset
    if nargin < 4
        INTERP = max(cell2mat(cellfun(@(x) size(x,2), D, 'UniformOutput', 0)));
    end

    % Normalize by length using interpolation
    nrm = @(d) [d ; 1 : length(d)]';
    L   = cellfun(@(x) getDim(interpolateOutline(nrm(x), INTERP), 1)', ...
        D, 'UniformOutput', 0);
    
else
    % Zero pad to longest length in dataset
    INTERP = max(cell2mat(cellfun(@(x) size(x,2), D, 'UniformOutput', 0)));
    padIt  = @(d) INTERP - size(d, 2);
    L      = cellfun(@(x) padarray(x', padIt(x), 0, 'pre')', ...
        D, 'UniformOutput', 0);
end

%% Handle width
if nrmW
    % Normalize widths to 1
    W = cellfun(@(x) x / max(x), L, 'UniformOutput', 0);
    W = fliplr(cat(1, W{:}));
    
else
    % Use raw widths
    W = fliplr(cat(1, L{:}));
    
end

%% Data regressors for Graham Schmidt
wids = cell2mat(cellfun(@(x) max(x), D, 'UniformOutput', 0));
lngs = cell2mat(cellfun(@(x) length(x), D, 'UniformOutput', 0));

end