function SUB = loadSubDirectories(DIN, dirName)
%% loadSubDirectories: Collect all sub-directories from root
% This is a helper script that creates a cell array of sub-directories from the
% inputted root directory DIN.
%
% Usage:
%   SUB = loadSubDirectories(DIN, dirName)
%
% Input:
%   DIN: root directory containing sub-directories
%   dirName: sub-directory from sub-directories containing image data
%
% Output:
%   SUB: cell array of sub-directories
%

%% Collect all sub-directories named dirName
% Get sub-directories from root directory
dins = dir2(DIN);

% Get all subdirectories named 'binary-masks'
% Doesn't work when sub-directories have multiple sub-directories
SUB = cell(1, numel(dins));
n = 1;
for din = dins'
    d = [din.folder '/' din.name];
    e = dir2(d);
    
    msks = e(cell2mat(arrayfun(@(x) strcmpi(x.name, dirName), ...
        e, 'UniformOutput', 0)));
    if ~isempty(msks)
        SUB{n} = [msks.folder '/' msks.name];
        n      = n + 1;
    end
end

% Remove empty cells
SUB = SUB(cell2mat(cellfun(@(x) ~isempty(x), SUB, 'UniformOutput', 0)));

end

