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
%         carrotStraightener(DIN, dirName, savData, savFigs, vis, par)
%
% Input:
%   DIN: path to root directory of image sub-directories
%   dirName: name of the sub-directory containing binary-masks
%   savData: boolean to save results in image directories, .mat, and .csv files
%   savFigs: boolean to save figures (vis parameter must be set to true)
%   vis: boolean to visualize output
%   par: boolean to run with parallel processing (1) or with normal loop (0)
%
% Output:
%   mlines: midlines 
%   cntrs: contours 
%   smsks: straightened masks 
%   pmsks: cell array of processed masks
%   tcrd: tip coordinates
%   dsts: sums of columns from straightened masks
%   fnames: filenames of images inputted to pipeline
%
% Example:
%   Run straightening pipeline on sub-directories named 'binary-masks' from
%   root directory, and save only output data:
%       % Root directory path and name of sub-directory of masks
%       din = '~/LabData/CarrotSweeper/z_datasets/masks_wi2019';
%       nm  = 'binary-masks';
%
%       % Run straightening pipeline from root directory with parellelization
%       [mlines, cntrs, smsks, pmsks, tcrds, dsts, fnames] = ...
%           carrotStraightener(din, nm, 1, 0, 1, 1);
%

%% Collect all sub-directories named dirName
% Get sub-directories from root directory
X = loadSubDirectories(DIN, dirName);

%% Normal run or with parallel processing
if par
    % Version that runs algorithm with parallel processing
    par = 0; % Can't run nested parfor loops
    parfor i = 1 : numel(X)
        x = X{i};
        [mlines{i}, cntrs{i}, smsks{i}, pmsks{i}, tcrds{i}, dsts{i}, ...
            fnames{i}] = carrotExtractor(x, vis, savData, savFigs, par);
    end
else
    % Run algorithm on all sub-directories and return data
    [mlines, cntrs, smsks, pmsks, tcrds, dsts, fnames] =  cellfun(@(x) ...
        carrotExtractor(x, vis, savData, savFigs, par), X, 'UniformOutput', 0);
end

end

