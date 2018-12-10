function fileNamesList=filesListFromInput(inputsList, flagGUI, filesFilter, setInptsSource, browserTitle)
%% filesListFromInput
% Returns cell array of file names, needed by various functions.
%
%% Syntax
%	fileNamesList=filesListFromInput(inputsList);
%   fileNamesList=filesListFromInput(inputsList, flagGUI)
%
%% Description
% This functions goal generate a cell array of file names, needed by the calling
%   functions. The input should be a a cell array with directories including the files or
%   files name. The function also suports input of a single directory name string, or a
%   single file name string. In the late two cases, each file is tested for existance.
%   Absolute file path is used thus replacing the relative path.
%   Alternatively, the user can choose the files or directories including fyles using the
%   OS explorer- by enabling the 'flagGUI' input.
%
%% Input arguments (defaults exist):
%	inputsList-    a path to the directory including files. A file name, or a cel array
%       of file names also supported.
%   flagGUI-      When enabled allowes the user to choose the files list using the
%       Explorer. Note files will be order accourding to their names (and not the order
%       you've clicked them).
%   filesFilter- a list (cell array) of extentions describing the file type user wishes to
%       detect. Files of types other then ones described in filesFilter will be ignored
%       and removed from the file names list.
%   setInptsSource- the input files source type- files list a directory or all files
%       located under the directory (including sub-folders). The user must choose one of
%       the three options {'Files list', 'Directory', 'Recursive Directory'}
%   browserTitle- a string used in the files/directories browser menu.
%
%% Output arguments
%	fileNamesList-    a cell array of file names, with absolute path.
%
%% Issues & Comments
%
%% Example I
%	fileNamesList=filesListFromInput(pwd);
%   fprintf('File names (+path) in current directory\n');
%   fprintf('%s\n', fileNamesList{:});
%
%% Example II
%   currDir=pwd;
%   cd( cat(2, matlabroot, '\toolbox\images\imdemos') );
%   fileNamesList=filesListFromInput( [], true );
%   cd(currDir);
%   fprintf( 'Chosen file names (+path) from Matlab images directory:\n' );
%   fprintf( '%s\n', fileNamesList{:} );
%
%% See also
%  - filesFullName
%
%% Revision history
% First version: Nikolay S. 2012-04-30.
% Last update:   Nikolay S. 2013-03-05.
%
% *List of Changes:*
% - 2013-03-05- The folderFiles function is used to get files recursivelly. Used by
%       setting 'nFolderDepth' input to 'Inf'.
% - 2012-11-04- An update that ignores files with extentions out of the filesFilter list.
% - 2012-08-28- Following inputs added: to reduce user clicking, when things (file types
%     and source) are well defined.
% - 2012-07-19- Empty input causes opening a files/directory explorer.
% - 2012-02-08- 'Unknown error occurred.' in fileattrib Matlab function taken care of.
%

%% Default params values
if nargin<2
    if nargin==0 % if no inputs supplied- force using explorer
        inputsList={};
    end
    if isempty(inputsList)
        flagGUI=true;
    else
        flagGUI=false; % by default, explorer will not be used
    end
end

if nargin < 3
    filesFilter=[];
end

if nargin < 4
    setInptsSource=[];
end

% when user did nor specify flag value, but specified further parameters- enable browser
if isempty(inputsList) && isempty(flagGUI)
    % if user has explicitly set flagGUI=false,  it will remain false,
    % (good for preventing unwanted user promts)
    flagGUI=true;
end

if isempty(flagGUI)
    flagGUI=false;
end

if exist('browserTitle', 'var')~=1
    browserTitle='Select input files';
end

if ischar(inputsList) % convert char to cell array
    inputsList={inputsList}; % consider using "cellstr" instead of {}
end

if (flagGUI)
    %% Allow the user to choose files/folder via explorer
    % Create a files filter aimed for user defined files, video files, image files, or all other files.
    if ~isempty(filesFilter)
        filesFilterFormated=sprintf('*.%s;',filesFilter{:});
        FilterSpec={filesFilterFormated,...
            cat(2, 'User defined files (', filesFilterFormated, ')');};
    else
        % get the file extentions of graphical and video formats supported by Matlab
        imageFormats=imformats;
        imageFormatsExtCell=cat(2, imageFormats.ext);
        imagesFilesFilter=sprintf('*.%s;',imageFormatsExtCell{:});
        videoFormats= VideoReader.getFileFormats();
        videosFilesExtList={videoFormats.Extension};
        videosFilesFilter=sprintf('*.%s;',videosFilesExtList{:});
        FilterSpec={ '*.*', 'All Files'; ...
            imagesFilesFilter, cat(2, 'Image Files (', imagesFilesFilter, ')');...
            videosFilesFilter, cat(2, 'Video Files (', videosFilesFilter, ')'); };
    end
    
    if isempty(inputsList)
        foldersList=inputsList;
    else
        % Devide inputsList to files and folders
        isFolder=cellfun(@isdir, inputsList);
        foldersList=inputsList(isFolder);
    end
    if isempty(foldersList)
        explStartDir=pwd; % start Explorer in current directory
    else
        explStartDir=folderFullPath( foldersList{end} ); % start Explorer in last of user folders
        inputsList{isFolder(end)}=[];
    end
    
    anotherDir='More';
    fileNames={};
    while ~strcmpi(anotherDir,'Finish')
        if exist( 'setInptsSource', 'var')==1 &&  any(strcmpi(setInptsSource,...
                {'Files list', 'Directory'}) ); 
            % user may choose one of the two supported browser types
            inptsSource=setInptsSource;
        else
            inptsSource = questdlg( 'Please choose files source:', browserTitle,...
                'Files list', 'Directory', 'Files list' );
        end
        
        switch(inptsSource)
            case{'Directory'}
                chosenDir = uigetdir(explStartDir, browserTitle);
                % store last opened directory, to start with it on next Explorer use.
                if isequal(chosenDir, 0) % If cancel was pressed
                    % skip iteration
                    continue;
                else % If a directory was chosen
                    explStartDir=chosenDir; % Use user selected folder
                end
                
                nFolderDepth = questdlg('Please choose folders depth:', browserTitle,...
                    '0-Current folder files', 'Custom value', 'Inf-All sub-folders files',...
                    'Inf-All sub-folders files');
                switch(nFolderDepth)
                    case('0-Current folder files')
                        nFolderDepth=0;
                    case('Custom value')
                        inputdlgPrompt='Enter folders depth value';
                        inputdlgTitle='Custom folders depth';
                        inputStr=inputdlg( inputdlgPrompt, inputdlgTitle, 1, {'1'} );
                        nFolderDepth=round( str2double(inputStr) );
                    otherwise
                        nFolderDepth=Inf;
                end     % switch(nFolderDepthUserChoise)

                fileNames=folderFiles(explStartDir, nFolderDepth, filesFilter);
            case{'Files list'}
                [fileName, pathName, ~] = uigetfile(FilterSpec, browserTitle,...
                    'MultiSelect', 'on', explStartDir);
                if iscell(fileName) || ischar(fileName) % cancel was not pressed
                    fileNames=cat(2, fileNames, strcat(pathName, fileName));
                    % store last opened directory, to start with it on next Explorer use.
                    explStartDir=pathName;
                end
        end % switch(inptsSource)
        
        anotherDir = questdlg({'Need to choose additional files?',...
            'Press ''More'', to choose additional files.',...
            'Press ''Finish'' to finish choosing inputs.'},...
            'Inputs files selection',...
            'More', 'Finish','Finish');
    end % while ~strcmpi(anotherDir,'Finish')
    
    % add GUI files, to users inputsList input
    if ~isequal( size(fileNames, 2), numel(fileNames) )
        fileNames=reshape(fileNames, 1, []);
    end
    
    inputsList=cat( 2, inputsList, fileNames ); 
    flagGUI=false; % this should be reconsidered
end % if (flagGUI)

% Devide inputsList to files and folders
inputsList=inputsList( ~cellfun(@isempty, inputsList) );
isFolder=cellfun(@isdir, inputsList);
foldersList=inputsList(isFolder); % Folders list

filesList=inputsList(~isFolder); % remove folders from the inputs list to get files list
% filesList=cellfun( @filesFullName, filesList, 'UniformOutput', false );
for iFile=1:length(filesList) % use full file path
    filesList{iFile}=filesFullName( filesList{iFile}, filesFilter, [], false );
end

isEmptyFile=cellfun( @isempty, filesList );
filesList=filesList( ~isEmptyFile );
if strcmpi( setInptsSource, 'Recursive Directory' )
    nFolderDepth=Inf;
else
    nFolderDepth=0;
end
folderFilesList=folderFiles(foldersList, nFolderDepth, filesFilter, flagGUI);

if ~isequal( size(folderFilesList, 2), numel(folderFilesList) )
    folderFilesList=reshape(folderFilesList, 1, []);
end     % if not(isequal( size(folderFilesList, 1), numel(folderFilesList) ))
if ~isequal( size(filesList, 2), numel(filesList) )
    filesList=reshape(filesList, 1, []);
end     % if not(isequal( size(folderFilesList, 1), numel(folderFilesList) ))

fileNamesList=cat(2, filesList, folderFilesList);

% % remore repeating files
% uniqueFilesList=unique(fileNamesList);
% if ~isequal( length(uniqueFilesList), length(fileNamesList) )
%     % If some file repetitons detected, remove obsolete elements
%     fileNamesList=uniqueFilesList;
% end