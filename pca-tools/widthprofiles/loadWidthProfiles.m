function FOUT = loadWidthProfiles(rootDir, maskDir, load_data)
%% loadWidthProfiles: load raw width profiles from directory path
% This function flips the widths left-to-right!
%
% TODO: 
% Write a class (perhaps Carrot?)  for these images that have methods to perform
% all the different  analyses from CarrotSweeper. This script will then create 
% instances of Carrot objects. 
%
% Usage:
%   [DSTS, IMGS, FNAMES] = loadWidthProfiles(rootDir, maskDir, load_data)
%
% Input:
%   rootDir: path to directory of root folder to straightened masks
%   maskDir: name of directory
%   load_data: only return path names [0] or compute images and profiles [1]
%
% Output:
%   FOUT: structure of outputs (only file paths if load_data set to 0)
%       paths: file paths to images
%       images: straightened masks used for getting width profiles
%       profiles: width profiles of straightened masks
%

%%
if nargin < 3
    load_data = 0; % Default to only return path names
end

%%
PATHS = loadSubDirectories(rootDir, maskDir);
STORE = imageDatastore(PATHS, 'IncludeSubfolders', 1, 'FileExtensions', '.png');

% Full file paths
FOUT = STORE.Files;

%% Load images and width profiles
% NOTE: very heavy in memory with large datasets!
if load_data
    FOUT.paths  = FOUT;
    FOUT.images = cellfun(@(x) x, STORE.readall, 'UniformOutput', 0);
    
    % Width profiles for all de-tipped images
    FOUT.profiles = cellfun(@(x) sum(logical(x)), IMGS, 'UniformOutput', 0);
    
end
end

