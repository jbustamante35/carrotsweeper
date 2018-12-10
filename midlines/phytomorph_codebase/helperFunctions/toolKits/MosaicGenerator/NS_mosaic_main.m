function outFile=NS_mosaic_main(varargin)
%
%% NS_mosaic_main
% Generates a mosaic of small tile images combined with a large main image.
%
%% Syntax
%  NS_mosaic_main; % Use previous run setting, or inputs via GUI
%  NS_mosaic_main('mainImgFile' ,mainImgFile, 'mosaicImgList', mosaicImgList,...
%     'sharpenEdges', sharpenEdges, 'unshAlpha', unshAlpha,...
%     'isGraySubs', isGraySubs, 'overlayRegionMat', overlayRegionMat, 'outFile', outFile);
%  any parameterrs can be ommited
%
%% Description
% This functions goal generate a mosaic image based on a single (usually large) main
%   image, and multiple small tile images. A large mosaic image will be randomly generated
%   from the multiple tile images. The resulting mosaic image will be combined with the
%   main image, forming the desired result. The ratio between images “strength” is defined
%   by the Overlay Ratio ranging between [0:1]. Smaller value will make the large main
%   image be dominant, while larger Overlay Ratio will make the mosaic image (formed from
%   tiles images) more obvious, at the expense of the main image. The user can (and
%   should) define regions with different Overlay Ratio values. For example when the  main
%   image is a portrait, the Overlay Ratio value around the face should be small (about
%   0.2), while it can be high (about 0.7) in other less important areas. The user can
%   also resize the final image, sharpen the small images, convert them to color (RGB) or
%   Gray images, if needed.
%
%% Input arguments (defaults exist):
% mainImgFile- the image to be used as mosaic base. Preferably a big, high resolution
%   image.
% mosaicImgList- a list of files and directories with images to be used as mosiac bulding
%   blocks- tiles. Small images will do.
% sharpenEdges- enhances small images, to make them salient.
% unshAlpha- Unsahrp alpha value [0:1]. Defines edges enhancement strength.
% isGraySubs- when enabled turn all small images to gray- so main images color will
%  remain unchanged.
% overlayRegionMat- the weight deteremening the influnce of sub images and main images in
%   the final mosiac [0:1]. Higher value will supress the main image, while lower value
%   will attenuate the main image.
% outFile- the resulting mosaic image file name and path. Default exists.
%
%% Output arguments
% outFile- the resulting mosaic image file name and path.
%
%% Issues & Comments
% The code is below my usual programing level, as features were added over time. It ought
% to be rewritten one day, but I must find the time. A classical example of "Spagetty Code"
%
%% Example
% NS_mosaic_main;
%
%% See also
%
%% Revision history
% First version: Nikolay S. 2011-01-15.
% Last update:   Nikolay S. 2014-08-14.
%
% *List of Changes:*+
% 2014-08-14:
%   - 'mosaicImages' updated by adding new files, removing obsolite ones, based on user previous and
%       current mosaic file choosing. So the same files are not uselsesly loaded. this can be a nice
%       way to generate a new mosiac form same mosiac Image Files, withotu actually reading them all
%       over again. just choose them as before, and feel free to change all the rest. The images are
%       loaded in the memory since the previous run!
% 2014-08-11:
%   - 'mosaicImages' resized to become Nx100, or 100xM, to save memory.
%   - mosaic images cropping prior to resizing updated. Cropping around center, preserving
%       acpect ratio.
% 2014-07-20:
%   - 'mosaicImages' varibale converted to persistent- to prevent reading the tile files over and
%   over again.
%   - an updates version of filesListFromInput is added, replacing older one.
% 2013-03-05:
%   - some documentation inprovements
%   - filesListFromInput changed to support recursive files listing via (slightly
%       modified) fuf.m
%       cortecy of Francesco di Pierro
% 2012-02-06:
%   - fileparts output variavles list changed to prevent error in new Matalb releases
%       (2011b and later).
%   - mosiac element choosing GUI modified- nor sqaure selection is supported, size is
%       presented on axis title.
%   - now user can choose directries and single files for the mosaic elements.
% 2012-02-14
%   - Blending of mosaic and main image is changed- now the overlay effects both.
%   - outFile naming convention and overwriting a file warning added.
%
%% Short Programs description
% Stages:
% 1)  Get main image, present it on screen.
% 2)  Suggest the user to resize it to different dimetions
% 3)  Ask user to mark mosiac image size with mouse on top of selected
%   image. Close the main image figure.
% 4)  Open file Explorer menu and ask the user to browse to a library with
%   the mosaic images.
% 5)  Open a GUI with following options:
%   5.1) RGB/Grayscale moisiac.
%   5.2) Mosiac elements aspect ratio (square or user defined).
%   5.3) Mosaic shapren strength
%   5.4) Mosaic and main image overlay ratio.
%       5.4.1) Let th user to define regions with different Overlay Ratios
% 6)  Mosiac image name (and file type).
% 7)  Read all mosiac images files/directories. Resize all images accourding to (3)
%   and (5.2). Calculate total resized elemnts area and compare it with
%   main image size. Define mosiac elemnts pool:
%   7.1) If image size bigger then mosiac elements area- randomly repeat some
%   elements- until desired numer of files mosiac elements will be choosen.
%   7.1) Otherwise- if image size smaller then mosiac elements area- randomly
%    choose some elements- until desired numer of files mosiac elements will
%    be choosen.
% 8)  Sharpen each pool element accourding to (5.3).
% 9)  Generate 2D/3D-(5.1)matrix of Inf, of size equal to main image.
% 10) Fill it with randomly chosen pool elemnts (randperm).
% 11) Search for remaining Infs and fill them with pool elements.
% 12) In case of 3D RGB in (5.1) repmat Grayscale matrix from (9) to RGB.
% 13) Combine mosaic and main image using (5.4).
% 14) Save and present resulting image.
% 15) Save mosiac generation setting under main image name with .mat file
%       extention.
%
% Notes- present waitbar with stage number description and progres bar of
% current stage.

%% TODO-
% save under main image mainFileName,
% better explain overlay,

%% load last used user parameters if exists, or default values otherwise
% mosaicImgList=pwd;
% currDir=pwd;
% addpath(currDir);
unshAlpha=0.2;
lastSettingFile='lastNSmosaic.mat';
% overlayRatio=0.5;

if nargin>0
    %% load command line user parameters
    for arg_ind=1:2:nargin
        eval([varargin{arg_ind},'=',varargin{arg_ind+1},';']);
    end
end
persistent mosaicImages;
persistent mosaicImgList;

% close all;clc;
imageFormats=imformats;
imageFormatsExtCell=cat(2, imageFormats.ext);

if exist(lastSettingFile, 'file')==2
    loadLastSetFile = menu('Program settings',...
        'Choose manually', 'Use last run params', 'Delete last run params');
else
    loadLastSetFile =1;
end

if (loadLastSetFile==2)
    %% Load data from previous run, with ony a few details to be changed by user
    load(lastSettingFile);
    
    hMainFig=figure;
    [mosaic_width, mosaic_heigth, main_img]=get_mosaic_size(mainImgFile);
    [main_img_rows, main_img_cols, ~]=size(main_img);
    needed_no_of_elements=ceil(main_img_rows/mosaic_heigth)*...
        ceil(main_img_cols/mosaic_width);
    close(hMainFig);
    
    % Basically- no need to reload- persistent variable is used !,
    % unless Matlab was closed, and the variable was deleted
    if exist('mosaicImages', 'var')==1 && isempty(mosaicImages)
        nFiles=length(mosaicImgList);
        mosaicImages=cell(1, nFiles);
        for iFile=1:nFiles
            currImg=imread( mosaicImgList{iFile} ); % read image file
            imgSize=size( currImg );    % resize, to x by 100 or 100 by x dims
            mosaicImages{iFile}=imresize(currImg, min( 100./imgSize([1, 2]) ));
            waitbarTimeRemaining2(iFile/nFiles, 'Loading & resizing mosiac files.');
        end
        
        mosaicProcImgs=prepare_mosiac_elements(mosaicImages,...
            mosaic_heigth, mosaic_width, needed_no_of_elements);
    end
    
    
else
    if (loadLastSetFile==3)
        delete(lastSettingFile);
        mosaicImages={};    % empty memory consuming persisten variable
        mosaicImgList={};   % Get read of irrelevant files list 
    end
    %% Or use various GUI menus to get inputs from user
    % Let the user choose file using browser- start with last used file
    if exist('mainImgFile', 'var')~=1
        mainImgFile=[];
    end
    mainImgFile=filesFullName(mainImgFile, imageFormatsExtCell, [],...
        'Select the Main Mosaic image file');
    assert( exist(mainImgFile, 'file')==2,...
        'Matlab can not acess this file (can be path/name problem).' );
    [mainFilePath, mainFileName, mainFileExt]= fileparts(mainImgFile);
    
    %% Or use various GUI menus to get inputs from user
    % Let the user choose directory using browser- start with last used file
    chosenMosaicFiles=filesListFromInput(mainFilePath, true, imageFormatsExtCell, [],...
        'Choose mosaic elements (files/folders)');

    if exist('mosaicImgList', 'var')~=1 || isempty(mosaicImgList)
        mosaicImgList=chosenMosaicFiles;
        % isEmptyFile=cellfun(@isempty, mosaicImgList, 'UniformOutput', false);
        % isEmptyFile=cat(1, isEmptyFile{:});
        % mosaicImgList=mosaicImgList(~isEmptyFile);
    else
        newFiles=setdiff(chosenMosaicFiles, mosaicImgList);
        [sharedFiles, isMosaic, ~]=intersect(mosaicImgList, chosenMosaicFiles, 'stable');
        mosaicImgList=cat(2, sharedFiles, newFiles);
        mosaicImages=mosaicImages(isMosaic); % remove obsolite elements
    end
    
    if exist('sharpenEdges', 'var')~=1
        choice = questdlg('Enable edges sharpening?', 'Sharpen selection:',...
            'Yes', 'No', 'Yes');
        if (strcmpi(choice, 'YES'))%(menu('Sharpen menu','On','Off')==1)
            sharpenEdges=true;
            inputdlg_opt.Interpreter='tex';
            unshAlpha=str2double(inputdlg({'Unsharp param \alpha val:'},...
                'Enter Unsharp params', 1, {'0.2'},inputdlg_opt));
        else
            sharpenEdges=false;
        end
    elseif sharpenEdges && exist('unshAlpha', 'var')~=1
        inputdlg_opt.Interpreter='tex';
        unshAlpha=str2double(inputdlg({'Unsharp param \alpha val:'},...
            'Enter Unsharp params', 1, {'0.2'},inputdlg_opt));
    end
    
    if exist('isGraySubs', 'var')~=1
        isGraySubs=false;
        choice = questdlg('Choose mosaic color', 'Color mode:',...
            'Gray', 'RGB', 'Gray');
        if (strcmpi(choice, 'GRAY')) % (menu('Mosiac color menu','Gray','RGB')==1)
            isGraySubs=true;
        end
    end
    
    %% Loading parameters is almost done, finish and start constructing mosaic image
    % present mosaic, resize if needed, let the user determine elemnt size
    hMainFig=figure;
    [mosaic_width, mosaic_heigth, main_img]=get_mosaic_size(mainImgFile);
    
    [main_img_rows, main_img_cols, ~]=size(main_img);
    needed_no_of_elements=ceil(main_img_rows/mosaic_heigth)*...
        ceil(main_img_cols/mosaic_width);
    
    % Read mosaic element images
    nFiles=length(mosaicImgList);
    if exist('newFiles', 'var')==1 && ~isempty(newFiles)
        nNewFiles=length(newFiles);
        iReadFile=nFiles+1-nNewFiles; 
        mosaicImages=cat( 2, mosaicImages, cell(1, nNewFiles) ); % add space for new images
    else
        iReadFile=1;
        mosaicImages=cell(1, nFiles);
    end
    for iFile=iReadFile:nFiles
        currImg=imread( mosaicImgList{iFile} ); % read image file
        imgSize=size( currImg );    % resize, to x by 100 or 100 by x dims
        mosaicImages{iFile}=imresize(currImg, min( 100./imgSize([1, 2]) ));
        waitbarTimeRemaining2( (iFile-iReadFile+1)/(nFiles-iReadFile+1),...
            'Loading & resizing new mosiac files.' );
    end
    
    %     % Prepare chosen images from mosiacing: crop and resize.
    %     mosaicProcImgs=prepare_mosiac_elements(mosaicImages, mosaic_heigth,...
    %         mosaic_width, needed_no_of_elements);
    
    
    overlayRegionChoise='Yes';
    inputdlg_opt.Interpreter='tex';
    overlayRegionMat=ones(main_img_rows,main_img_cols, 'double');
    overlayRatio=str2double(...
        inputdlg({'Overlay ratio of non selected regions [0:1]'},...
        'Combine main and mosaic images menu', 1, {'0.5'}, inputdlg_opt));
    if ~isempty(overlayRatio)
        overlayRatio=min(overlayRatio, 1); % Overlay ratio must be under 1
        overlayRatio=max(0, overlayRatio); % Overlay ratio must be over 0
        overlayRegionMat=overlayRatio*overlayRegionMat;
    end
    
    %% Defining regions with different Overlay values
    hold on;
    while ~strcmpi(overlayRegionChoise, 'No')
        overlayRegionChoise=questdlg('Define another region with diferent Overlay?',...
            'Overlay regions selection:', 'Yes', 'No', 'Yes');
        
        if strcmpi(overlayRegionChoise, 'No')
            continue;
        end
        title('Select by left clicking mouse and marking the Region.',...
            'FontSize',18, 'Color', [1,0,0]);
        imfreehandH=imfreehand;
        selMask=createMask(imfreehandH);
        posRegion = getPosition(imfreehandH);
        posRegion=cat(1, posRegion, posRegion(1,:));
        delete(imfreehandH);
        
        plot(posRegion(:,1), posRegion(:,2), 'Color', rand(1,3), 'LineWidth', 2 );
        overlayRatio=str2double( inputdlg({'Overlay ratio [0:1]'},...
            'Combine main and mosaic images menu', 1, {'0.5'},...
            inputdlg_opt) );
        if ~isempty(overlayRatio) % if "cancel was not chosen"
            overlayRatio=min(overlayRatio, 1); % Overlay ratio must be under 1
            overlayRatio=max(0, overlayRatio); % Overlay ratio must be over 0
            
            overlayRegionMat(selMask)=overlayRatio;
        end
        title('Main image, subject to mosaicing', 'FontSize', 14,...
            'Color', [0,0,0]);
    end % while ~strcmpi(overlayRegionChoise, 'FINISH.')
    hold off;
    close(hMainFig);
    
    %% Smooth the sharp transitions betwen the user defined regions
    smoothFilt=ones(ceil(mosaic_heigth/2), ceil(mosaic_width/2));
    smoothFilt=smoothFilt/sum(smoothFilt(:));
    overlayRegionMat=filter2(smoothFilt, overlayRegionMat, 'same');
    
end  % if(loadPrevSetFile==1)

%% Prevent overwriting exisitng files
if exist('outFile', 'var')~=1
    outFile=[mainFilePath, filesep, 'mosaicNS_' mainFileName, mainFileExt];
end % if exist('outFile', 'var')~=1

if exist(outFile, 'file')==2 % Promtpt user not to override existing file
    choice = questdlg('File exist, overwrite?', 'Warning.',...
        'Yes', 'No', 'Yes');
    if strcmpi(choice, 'NO')
        [fileOutPath, ~, fileOutExt]=fileparts(outFile);
        [~, mainFileName, ~]= fileparts(mainImgFile);
        outFile=[fileOutPath, filesep, 'mosaicNS_' mainFileName,'_',...
            datestr(now,'yyyymmmdd_HH-MM-SS'), fileOutExt];
    end
end % if exist(outFile, 'file')==2 % Prompt user not to override existing file

dlg_cell=inputdlg({'Save result to file:'}, 'Enter file name', 1, {outFile});
outFile=dlg_cell{1};

% save current run settings
save(lastSettingFile, 'mainImgFile', 'outFile', 'unshAlpha',...
    'sharpenEdges', 'isGraySubs', 'overlayRegionMat', 'mosaicImgList');

mosaicProcImgs=prepare_mosiac_elements(mosaicImages, mosaic_heigth,...
    mosaic_width, needed_no_of_elements);
no_of_elems=length(mosaicProcImgs);
% Sharpen via Unsharp+
if sharpenEdges
    h = waitbar(0,'Sharpening mosaic images.');
    
    for mosaic_images_ind=1:no_of_elems
        mosaicProcImgs{mosaic_images_ind}=imfilter(mosaicProcImgs{mosaic_images_ind},...
            fspecial('unsharp',unshAlpha),'replicate');
        waitbar(mosaic_images_ind/no_of_elems,h);
    end
    close(h);
end

% Randiomly sort mosiac elements
mosaic_el_perm=permute_mosaic_elements(needed_no_of_elements,no_of_elems);

% define final image mosiac part
res_mosaic=zeros(ceil(main_img_rows/mosaic_heigth)*mosaic_heigth,...
    ceil(main_img_cols/mosaic_width)*mosaic_width,3,'uint8');

% locate the mosiac elemnts properly
h = waitbar(0, 'Forming mosaic.');
mosaic_elem_ind=1;
for rows_ind=1:size(res_mosaic,1)/mosaic_heigth
    for cols_ind=1:size(res_mosaic,2)/mosaic_width
        curr_mosaic_element=mosaicProcImgs{mosaic_el_perm(mosaic_elem_ind)};
        
        res_mosaic( (1+(rows_ind-1)*mosaic_heigth):rows_ind*mosaic_heigth,...
            (1+(cols_ind-1)*mosaic_width):cols_ind*mosaic_width,:)=curr_mosaic_element;
        mosaic_elem_ind=mosaic_elem_ind+1;
    end
    waitbar(rows_ind/(size(res_mosaic,1)/mosaic_heigth),h);
end
close(h);

% Crop to original image size
res_mosaic=res_mosaic(1:main_img_rows, 1:main_img_cols, :);

% convert to Gray
if isGraySubs
    res_mosaic=repmat(rgb2gray(res_mosaic),[1,1,3]);
end

overlayRegionMat=repmat(overlayRegionMat,[1,1,3]);
if size(main_img,3)==1
    main_img=repmat(main_img, [1,1,3]);
end
if ~isequal( size(overlayRegionMat), size(main_img) )
    overlayRegionMat=imresize(overlayRegionMat, [size(main_img, 1), size(main_img, 2)],...
        'nearest' );
end
% combine original image with mosiac, to get the mosaic final image.
main2mosaicG=1; % values under 1 will preserve strength of main image
final_image=uint8( (1-overlayRegionMat*main2mosaicG).*double(main_img)+...
    overlayRegionMat.*double(res_mosaic) );

% present final image
figure;
imshow(final_image);

% save final image
imwrite(final_image, outFile);