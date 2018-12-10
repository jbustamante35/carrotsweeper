function [] = main(varargin)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % define inpath
    if nargin == 0
        inFilePath = uigetdir();
        inFilePath = [inFilePath filesep];
    else
        inFilePath = varargin{1};
    end
    % define outpath
    oPath = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/massExtraction/';
    mkdir(oPath)
    % define flag path
    fPath = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/massExtraction/flagPath/';
    mkdir(fPath);
    % define errorpath
    ePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/massExtraction/errorPath/';
    mkdir(ePath);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% scan for new images
    FileList = {};
    FileExt = {'tiff','TIF','tif'};
    verbose = 1;
    SET = sdig(inFilePath,FileList,FileExt,verbose);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% remove non 61
    ridx = [];
    for e = 1:numel(SET)
        if numel(SET{e}) ~= 61
            ridx(e) = 1;
        else
            ridx(e) = 0;
        end
    end
    SET(find(ridx)) = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:numel(SET)
        tic
        try
            flagFile = [fPath strrep(fileparts(SET{e}{1}),'/','--') '.csv'];
            if ~exist(flagFile)
                isolateKernels(SET{e},1,oPath);
                %[data{e}] = sampleContours(SET{e},0,oPath);
                csvwrite(flagFile,1);
            else
                fprintf(['Skipping:' SET{e}{1} '\n']);
            end
            
        catch ME
            [pth nm ext] = fileparts((ME.stack(1).file));
            errorFile = [nm '~' num2str(ME.stack(1).line)];
            errorFile = [ePath errorFile '~' '{' ME.message '}~' strrep(fileparts(SET{e}{1}),'/','--') '.csv'];
            csvwrite(errorFile,1);            
        end
        
        
        %isolateKernels(SET{e},1,[]);
        toc
    end
end
%{
    main('/mnt/spaldingdata/Takeshi/allMaizeMovies/');
%}