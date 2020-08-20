function [DSTS, IMGS, FPATHS] = loadWidthProfiles(rootDir, maskDir)
%% loadWidthProfiles: load raw width profiles from directory path
% This function flips the widths left-to-right!
%
% Usage:
%   [DSTS, IMGS, FNAMES] = loadWidthProfiles(rootDir, maskDir)
%
% Input:
%   rootDir: path to directory of root folder to straightened masks
%   maskDir: name of directory
%
% Output:
%   DSTS: width profiles of straightened masks
%   IMGS: straightened masks used for getting width profiles
%   FNAMES: filenames of images
%

%%
if nargin < 1
    IMGS = '/home/jbustamante/Dropbox/EdgarSpalding/labdata/rawdata/scott_brainard/images';
    SET  = 'gwas-curated';
    
    rootDir = sprintf('%s/%s', IMGS, SET);
    maskDir = 'straight-masks';
end

%%
PATHS = loadSubDirectories(rootDir, maskDir);
STORE = imageDatastore(PATHS, 'IncludeSubfolders', 1, 'FileExtensions', '.png');
% IMGS  = cellfun(@(x) x, STORE.readall, 'UniformOutput', 0);
IMGS = [];

% Width profiles for all de-tipped images
% DSTS = cellfun(@(x) sum(logical(x)), IMGS, 'UniformOutput', 0);
DSTS = [];

% Full file paths
FPATHS = STORE.Files;

end

