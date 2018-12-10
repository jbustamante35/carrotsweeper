function [FileList] = removeProcessedFiles(FileList,user,tissueType,plantType)
    if nargin == 3
        plantType = 'maize';
    end
    % glue the plantType to the tissueType
    analysisType = [plantType '_' tissueType];
    % check to see if data is already analyzed
    switch analysisType
        case 'maize_ears'
            %resultsPath = ['/mnt/spaldingdata/nate/mirror_images/maizeData/' user '/return/earData/output/'];
            resultsPath = ['/home/nate/iplant/return/earData/'];
            sfield = '_result';
        case 'maize_cobs'
            %resultsPath = ['/mnt/spaldingdata/nate/mirror_images/maizeData/' user '/return/cobData/output/'];
            resultsPath = ['/home/nate/iplant/return/cobData/'];
            sfield = '_result';
        case 'maize_kernels'
            %resultsPath = ['/mnt/spaldingdata/nate/mirror_images/maizeData/' user '/return/kernelData/output/'];
            resultsPath = ['/home/nate/iplant/return/kernelData/'];
            sfield = '';
        case 'maize_seedlings'
            resultsPath = ['/home/nate/iplant/return/seedlingData/'];
            sfield = '{PlantNumber_All}{Phenotype_Image}';
        case 'carrot_wholes';
            resultsPath = ['mnt/spaldingdata/nate/mirror_images/carrotData/' user '/return/wholeData/output/'];
            sfield = '_result';
        case 'arabidopsis_rootStraight'
            resultsPath = ['mnt/spaldingdata/nate/mirror_images/arabidopsisData/' user '/return/rootStraightKinematics'];
            sfield = '_result';
    end
    
    % check to see if resulting image is present
    rm = [];
    for e = 1:numel(FileList)
        [pth nm ext] = fileparts(FileList{e});
        ext = '.tif';
        v1 = keyLookup(FileList{e},'Plot');
        if strcmp(v1,'126')
            nm
        end
        oFile = [resultsPath nm sfield ext];
        %rm(e) = exist(oFile);
    end
    
    % remove the already processed files
    rm = logical(rm);
    FileList(rm) = [];
end