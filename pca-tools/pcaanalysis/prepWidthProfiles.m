function [W, wids, lngs] = prepWidthProfiles(D, mth, I)
%% prepWidthProfiles: prepare width profiles for PCA
%
%
% Usage:
%   [W, wids, lngs] = prepWidthProfiles(D, mth)
%
% Input:
%   D: cell array of width profiles
%   mth: method of prep to use [zeropad|normalize]
%   I: size to interpolate length to
%
% Output:
%   W: width profiles stacked using requested method
%

if nargin < 3
    I = max(cell2mat(cellfun(@(x) size(x,2 ), D, 'UniformOutput', 0)));
end

switch mth
    case 'zeropad'
        padIt   = @(d) I - size(d, 2);
        W      = cellfun(@(x) padarray(x', padIt(x), 0, 'pre')', ...
            D, 'UniformOutput', 0);
        W      = fliplr(cat(1, W{:}));
        
    case 'normalize'
        % Interpolate by longest length
        nrm = @(d) [d ; 1 : length(d)]';
        NL  = cellfun(@(x) getDim(interpolateOutline(nrm(x), I), 1)', ...
            D, 'UniformOutput', 0);
        
        % Normalize by widths
        NW = cellfun(@(x) x / max(x), NL, 'UniformOutput', 0);                       
        W  = fliplr(cat(1, NW{:}));
        
        % Data regressors
        wids = cell2mat(cellfun(@(x) max(x), D, 'UniformOutput', 0));
        lngs = cell2mat(cellfun(@(x) length(x), D, 'UniformOutput', 0));
        
    case 'tips'
        
        
    otherwise
        fprintf(2, 'Method must be [zeropad|normalize]\n');
        W = [];
end

end