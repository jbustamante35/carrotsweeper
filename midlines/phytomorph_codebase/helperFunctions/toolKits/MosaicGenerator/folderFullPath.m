function fullPathFolder=folderFullPath(folderString)
%% folderFullPath
% The function attempts to find a full path of a partial folder path, filling the missing
%   gaps.
%
%% Syntax
%  fullPathFolder=folderFullPath(folderString);
%
%% Description
% The function uses Matlab build in commands "fileattrib" and "what" to get the file
%   details, one of which is the files full path we desire.
%
%% Input arguments (defaults exist):
%  inFile- input file name. inFile must include file name. File path may be ommited, if
%     file is in Matlab path. File extention can be mmited for simplisuty or laziness.
%  filesExtList- a list of extentions describing the file type we wish to detect. Default
%     file types are Graphical or Videos.
%  dlgTitle- a string used in the files explorer menu.
%  isErrMode- a logical variable defining fnction behaviour in case of non existent file.
%     When enabled- an error messge will be issued (default behaviour). When disabled, an
%     empty name will be returned, without an error
%
%% Output arguments
%   fullFileName-  a full file name (path+file name+extention).
%
%% Issues & Comments
% "fileattrib" command fails sometimes for an unknown reason, therefore slower "what"
%   command is used
%
%% Example
%
%% See also
% - folderSubFolders
%
%% Revision history
% First version: Nikolay S. 2012-05-01.
% Last update:   Nikolay S. 2012-11-14.
%
% *List of Changes:*

%% prepare folderString folder for further processing
% get full folder path, in case relative/partial path was used
fullPathFolder=folderString;

if isdir(folderString)
    [stats, currFolderAttr]=fileattrib(folderString);
    % sometimes "fileattrib" fails witout any explanation
    if ( stats && ~strcmpi(currFolderAttr, 'Unknown error occurred.') )
        fullPathFolder=currFolderAttr.Name;
    else
        % if folder exists but "fileattrib" function failed use slower "what"
        currFolderAttr=what(folderString);
        fullPathFolder=currFolderAttr(end).path; % index can also be 1
    end	% if ( stats && ~strcmpi(currFolderAttr, 'Unknown error occurred.') )
end	% if isdir(folderString)
