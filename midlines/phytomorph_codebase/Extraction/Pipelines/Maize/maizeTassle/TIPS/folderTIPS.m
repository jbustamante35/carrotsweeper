function [  ] = folderTIPS( directory )
size(directory)
for e = 1:numel(directory)
    directory{e}
end
%UNTITLED Runs TIPS on a set of images in a directory
%   Will look for images in the provided directory.  Only pairs of images
%   named <imgName>.jpg and <imgName>_background.jpg will be processed.

%{
directory = ['./' directory{1}]
outPath = './output/';


fileList = dir(fullfile(directory, '/*.jpg'));
fileList = {fileList.name}.';
%}
outPath = './output/';
mkdir(outPath);
fileList = directory';
[tp,tn,text] = fileparts(fileList{1});
tp = ['./' tp filesep]
fidx = strfind(tp,filesep);
directory = tp(1:fidx(2));
directory
size(fileList)

for e = 1:numel(fileList)
    [p,fileList{e},ext] = fileparts(fileList{e});
    fileList{e} = [fileList{e} ext]
    %CMD = ['mkdir -p ' outPath p filesep];
    %fprintf(['command: ' CMD '\n']);
    %system(CMD);
end

findBackground = @(str) cellfun(@(c) ~isempty(c), regexp(fileList, str, 'once'));
findForeground = @(str) cellfun(@(c) isempty(c), regexp(fileList, str, 'once'));
getBase = @(str, sub) strrep(str, sub, '');

% Get background and foreground names from file list
backgrounds = fileList(findBackground('.+_background\.jpg'));
foregrounds = fileList(findForeground('.+_background\.jpg'));

% Make cell arrays of matching foreground/background images
foreBase = cellfun(@(c) getBase(c,'.jpg'), foregrounds, 'UniformOutput', false);
backBase = cellfun(@(c) getBase(c, '_background.jpg'), backgrounds, 'UniformOutput', false);

common = intersect(foreBase, backBase);
foregroundsToKeep = strcat(common, '.jpg');
backgroundsToKeep = strcat(common, '_background.jpg');
baseToKeep = common;

% Get backgrounds w/o a foreground and vice versa
backNoFore = ~ismember(backBase, common);
backgroundsToSkip = backgrounds(backNoFore);

foreNoBack = ~ismember(foreBase, common);
foregroundsToSkip = foregrounds(foreNoBack);

% Print warning about backgrounds w/o foregrounds and vice versa
if ~isempty(foregroundsToSkip)
    fprintf('NOTICE: The following foregrounds have no background and will be skipped:\n\n');
    for i = 1:size(foregroundsToSkip, 1)
        fprintf(foregroundsToSkip{i});
        fprintf('\n');
    end
    fprintf('\n');
end

if ~isempty(backgroundsToSkip)
    fprintf('NOTICE: The following backgrounds have no foreground and will be skipped:\n\n');
    for i = 1:size(backgroundsToSkip, 1)
        fprintf(backgroundsToSkip{i});
        fprintf('\n');
    end
    fprintf('\n');
end

foregroundsToKeep
fprintf(['start running on:' num2str(size(foregroundsToKeep, 1)) ':images\n']);
% Run TIPS on matched foreground/background images
for i = 1:size(foregroundsToKeep, 1)
    foreground = fullfile(directory, foregroundsToKeep{i});
    background = fullfile(directory, backgroundsToKeep{i});
    out = fullfile(outPath, baseToKeep{i});
    class(foreground)
    foreground
    TIPS(foreground, background, out);
end
fprintf(['end running on:' num2str(size(foregroundsToKeep, 1)) ':images\n']);
close all

