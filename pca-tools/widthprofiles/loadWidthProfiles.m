function FOUT = loadWidthProfiles(rootDir, maskDir, load_data, img)
%% loadWidthProfiles: load raw width profiles from directory path
% This function flips the widths left-to-right!
%
% TODO:
% Write a class (perhaps Carrot?)  for these images that have methods to perform
% all the different  analyses from CarrotSweeper. This script will then create
% instances of Carrot objects.
%
% Usage:
%   FOUT = loadWidthProfiles(rootDir, maskDir, load_data, img)
%
% Input:
%   rootDir: path to directory of root folder to straightened masks
%   maskDir: name of directory
%   load_data: return path names [0], images [1], or width profiles [2]
%   img: image data type (default '.png')
%
% Output:
%   FOUT: structure of outputs (only file paths if load_data set > 0)
%       paths: file paths to images
%       images: straightened masks used for getting width profiles
%       profiles: width profiles of straightened masks
%

%% Parse inputs
switch nargin
    case 2
        load_data = 0;      % Default to only return path names
        img       = '.png'; % Default to search for png images
    case 3
        img       = '.png'; % Default to search for png images
    case 4
    otherwise
        fprintf(2, 'Error with inputs (%d)\n', nargin);
        FOUT = [];
        return;
end

%%
PATHS = loadSubDirectories(rootDir, maskDir);
STORE = imageDatastore(PATHS, 'IncludeSubfolders', 1, 'FileExtensions', img);

%% Load images and width profiles
% NOTE: very heavy in memory with large datasets!
if load_data
    FOUT.paths = STORE.Files;
    IMGS       = cellfun(@(x) x, STORE.readall, 'UniformOutput', 0);
    
    switch load_data
        case 1
            % Return just the straightened images
            FOUT.images = IMGS;
        case 2
            % Return the width profiles
            FOUT.profiles = cellfun(@(x) sum(logical(x)), ...
                IMGS, 'UniformOutput', 0);
        case 3
            % Return both images and profiles [heavy on RAM!]
            FOUT.images   = IMGS;
            FOUT.profiles = cellfun(@(x) sum(logical(x)), ...
                IMGS, 'UniformOutput', 0);
        otherwise
            fprintf(2, 'Incorrect input %d [1|2]\n', load_data);
    end
else
    % Just return filepaths
    FOUT = STORE.Files;
end
end

