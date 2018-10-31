function [out] = isolateRoots_overStack(stack,outPath,outName,toWrite,NP,NPK,SNIP,disp,numToProcess)
    % stack - the file names of the images    
    % outPath - the location to write the results to
    % outName - the name of the files to save the results as
    % toWrite - flag for writing csv files    
    % NP - number of points for tip angle measurement
    % NPK - number of points to measure kurvature over
    % SNIP - number of points to SNIP from midline
    % disp - to display
    % numImagesToProcess - if less than zero then process all
    
    if isdeployed
        fprintf(['Converting inputs \n']);
        % convert inputs
        stack = cellStringConvert(stack);
        toWrite = str2num(toWrite);
        NP = str2num(NP);
        NPK = str2num(NPK);
        SNIP = str2num(SNIP);
        disp = str2num(disp);
        numToProcess = str2num(numToProcess);
    end
    
    %%%%%%%%%%%%%%%%%%
    % loop over stack for extraction of midline
    if numToProcess < 0
        numToProcess = numel(stack);
    end
    
    
    %%%%%%%%%%%%%%%%%%
    % sor the images
    for e = 1:numel(stack)
        [p,n,ext] = fileparts(stack{e});
        NUM(e) = str2num(n);
    end
    [~,sidx] = sort(NUM);
    stack = stack(sidx);
    
    parfor e = 1:numToProcess
        tic
        fprintf(['Procesing filename:' stack{e} '\n'])
        out{e} = isolateRoots(stack{e},disp,outPath,num2str(e));
        fprintf(['Done with :' num2str(e) ':' num2str(numel(stack)) ':' num2str(toc) '\n']);
    end
    
    %%%%%%%%%%%%%%%%%%
    % save to mat file
    [pth,nm,ext] = fileparts(stack{1});
    %fidx = strfind(pth,filesep);    
    %pth = pth(fidx(end-DEPTH)+1:end);
    %pth = strrep(pth,filesep,'----');
    %matOutPath = [outPath 'mat/'];
    matOutPath = [outPath];
    mkdir(matOutPath);
    fileName = [matOutPath outName '.mat'];
    %%%%%%%%%%%%%%%%%%%%%
    % generate output name
    fileName = strrep(fileName,'SLASH','--');
    fileName = strrep(fileName,'SPACE',' ');
    save(fileName,'out');
    
    
    if toWrite
        % write out data
        %csvOutPath = [outPath 'csv/'];
        csvOutPath = [outPath];
        mkdir(csvOutPath);
        writeData(fileName,[csvOutPath outName],NP,SNIP,NPK);
    end
    close all
end

%{
%%%%%%%%%%%%%%%%%%%%%%%%
% run single instance of this program
inPath = '/mnt/spaldingdata/steve/Gravitropism/5-11-15/eir1;b4_cam1/';
inPath = '/mnt/spaldingdata/steve/Gravitropism/5-18/eir1-5_cam4/';
inPath = '/mnt/spaldingdta/nate/mirror_images/otherLab_gravitropism_mid/shih-heng/Gaz8/';
inPath = '//home/nate/Downloads/TR11A/';
%%%%%%%%%%%%%%%%%%%%%%%%
% set file ext
FileExt = {'tiff','TIF','tif','JPG'};         % the file extensions to scan for    
% set verbose
verbose = 1;
% perform scan
fileSet = gdig(inPath,{},FileExt,verbose);
toWrite = 1;
NP = 20;
NPK = 20;
SNIP = 20;
disp = 1;
numToProcess = -1;
[out] = isolateRoots_overStack(fileSet,'/mnt/snapper/nate/return/','forGrant',1,NP,NPK,SNIP,disp,numToProcess);

%}