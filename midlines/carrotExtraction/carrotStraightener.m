function [mlines, cntrs, smsks, pmsks, tcrds, dsts, fnames] = carrotStraightener(DIN, dirName, savData, savFigs, vis, par)
%% carrotStraightener: function to run carrotExtractor on sub-directories
% This is a script that simply takes in a root directory and identifies all the
% sub-directories matching the dirName parameter and runs the carrotExtractor
% algorithm to run the straightening pipeline.
%
% [Insert more description of the pipeline here]
%
% [Insert information about what data is outputted here]
%
% Usage:
%   [mlines, cntrs, smsks, pmsks, tcrds, dsts, fnames] = ...
%         carrotStraightener(DIN, dirName, savData, savFigs, vis)
%
% Input:
%   DIN: path to root directory of data to analyze
%   dirName: name of the sub-directory containing binary-masks
%   savData: boolean to save outputted contour, midline, straightened mask data
%   savFigs: boolean to save figures if the vis parameter is set to true
%   vis: boolean to show mask with midline/contour and straightened mask
%
% Output:
%   mlines: cell array of midlines for each sub-directory of images
%   cntrs: cell array of contours for each sub-directory of images
%   smsks: cell array of straightened masks for each sub-directory of images
%   pmsks: cell array of processed masks for each sub-directory of images
%   tcrd: cell array of tip coordinates
%   dsts: cell array of distance transform values along midline
%   fnames: cell array of filenames of images
%
% Example:
%   Run straightening pipeline on sub-directories named 'binary-masks' from
%   root directory, and save only output data:
%       % Root directory path and name of sub-directory of masks
%       din = '~/LabData/CarrotSweeper/z_datasets/masks_wi2019';
%       nm  = 'binary-masks';
%
%       % Run straightening pipeline from root directory
%       [mlines, cntrs, smsks, pmsks, tcrds, dsts] = ...
%           carrotStraightener(din, nm, 1, 0, 1);
%

%% Collect all sub-directories named dirName
% Get sub-directories from root directory
X = loadSubDirectories(DIN, dirName);

%% Run algorithm on all sub-directories and return data
% [mlines, cntrs, smsks, pmsks, tcrds, dsts, fnames] =  cellfun(@(x) ...
%     carrotExtractor(x, vis, savData, savFigs), X, 'UniformOutput', 0);

%% Normal run or with parallel processing
if nargin > 5
    if par
        % Version that runs algorithm with parallel processing
        parfor i = 1 : numel(X)
            x = X{i};
            [mlines{i}, cntrs{i}, smsks{i}, pmsks{i}, tcrds{i}, dsts{i}, ...
                fnames{i}] = carrotExtractor(x, vis, savData, savFigs);
        end
    else
        % Run algorithm on all sub-directories and return data
        [mlines, cntrs, smsks, pmsks, tcrds, dsts, fnames] =  cellfun(@(x) ...
            carrotExtractor(x, vis, savData, savFigs), X, 'UniformOutput', 0);
    end
end

end

