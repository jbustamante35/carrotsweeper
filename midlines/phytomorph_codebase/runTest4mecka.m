function [] = runTest4mecka()
    % warnings will not be printed on screen.
    warning off
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                runTest4mecka.m is to download and run test analysis.(Inputs are relative to 1200dpi)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Prepare test run for mecka
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        fprintf(['Runing mecka test: Designed on May 31 2016. \n']);
        fprintf(['Test will download and run on ear, cob and kernel data and run analysis \n']);
        [rootPath,nm,ext] = fileparts(which('runTest4mecka'));
        fprintf(['root test location is created at ' rootPath ' \n']);
        %%%%%%%%%%%%%%%%%%%%%%%
        % make input directory
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['create testImage directory \n']);
        inputdir = [rootPath, filesep, 'Tests', filesep, 'testImages'];
        CMD = ['mkdir -p ' inputdir];
        system(CMD);
        fprintf(['testImage directory is created at ' inputdir ' \n']);
        %%%%%%%%%%%%%%%%%%%%%%%
        % make output directory
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['create testResult directory \n']);
        outputdir = [rootPath, filesep, 'Tests', filesep, 'testResults'];
        CMD = ['mkdir -p ' inputdir];
        system(CMD);
        fprintf(['testResult directory is created at ' outputdir ' \n']);
        %%%%%%%%%%%%%%%%%%%%%%%
        % download test images
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['********************************************\n']);
        fprintf(['Starting download of test images from iPlant\n']);
        fprintf(['********************************************\n']);
        earImage = [inputdir, filesep, 'earImage.tif'];
        cobImage = [inputdir, filesep, 'cobImage.tif'];
        kernelImage = [inputdir, filesep, 'kernelImage.tif'];
        fprintf(['start download earImage at ' earImage ' \n']);
        urlwrite('http://davos.cyverse.org/irods-rest/rest/fileContents/iplant/home/nmiller/publicData/earImage.tif?ticket=nwsqWo4WIHpdB0B', earImage);
        fprintf(['end download earImage \n']);
        fprintf(['start download cobImage at ' cobImage ' \n']);
        urlwrite('http://davos.cyverse.org/irods-rest/rest/fileContents//iplant/home/nmiller/publicData/cobImage.tif?ticket=l1dt2blAYjNqHHd', cobImage);
        fprintf(['end download cobImage \n']);
        fprintf(['start download kernelImage at ' kernelImage ' \n']);
        urlwrite('http://davos.cyverse.org/irods-rest/rest/fileContents//iplant/home/nmiller/publicData/kernelImage.tif?ticket=nwgnXuIM1L7rpLM', kernelImage);
        fprintf(['end download kernelImage \n']);
        fprintf(['********************************************\n']);
        fprintf(['Ending download of test images from iPlant\n']);
        fprintf(['********************************************\n']);
        %%%%%%%%%%%%%%%%%%%%%%%
        % set path for runTest4mecka
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf('set path for test run \n');
        % add path: including subdirectories
        addpath(genpath(rootPath));
        % remove path: excluding .git
        gitPath = [rootPath, filesep, '.git'];
        rmpath(genpath(gitPath));
        fprintf([rootPath, ' and subdirectories added as path for test run \n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Run test analysis for mecka
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        fprintf(['********************************************\n']);
        fprintf(['Starting analysis of test images from iPlant\n']);
        fprintf(['********************************************\n']);
        %%%%%%%%%%%%%%%%%%%%%%%
        % run test analysis for ear
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['****************\n']);
        fprintf(['Starting analysis of ear image.\n']);
        fprintf(['****************\n']);
        fprintf(['create ear directory \n']);
        oPathE = [outputdir, filesep, 'Ear'];
        mkdir(oPathE);
        fprintf(['ear directory is created at ' oPathE ' \n']);
        mecka('e',earImage,3,oPathE,'remotePath',1,1,1200,1);
        fprintf(['****************\n']);
        fprintf(['Ending analysis of ear image.\n']);
        fprintf(['****************\n']);
        %%%%%%%%%%%%%%%%%%%%%%%
        % run test analysis for cob
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['****************\n']);
        fprintf(['Starting analysis of cob image.\n']);
        fprintf(['****************\n']);
        fprintf(['create cob directory \n']);
        oPathC = [outputdir, filesep, 'Cob'];
        mkdir(oPathC);
        fprintf(['cob directory is created at ' oPathC ' \n']);
        mecka('c',cobImage,3,oPathC,'remotePath',1,1,1200,1);
        fprintf(['****************\n']);
        fprintf(['Ending analysis of cob image.\n']);
        fprintf(['****************\n']);
        %%%%%%%%%%%%%%%%%%%%%%%
        % run test analysis for kernel
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['****************\n']);
        fprintf(['Starting analysis of kernel image.\n']);
        fprintf(['****************\n']);
        fprintf(['create kernel directory \n']);
        oPathK = [outputdir, filesep, 'Kernel'];
        mkdir(oPathK);
        fprintf(['kernel directory is created at ' oPathK ' \n']);
        mecka('k',kernelImage,3,oPathK,'remotePath',1,1,1200,1);
        fprintf(['****************\n']);
        fprintf(['Ending analysis of kernel image.\n']);
        fprintf(['****************\n']);
        fprintf(['********************************************\n']);
        fprintf(['Ending analysis of test images from iPlant\n']);
        fprintf(['********************************************\n']);
    catch ME
        close all;
        getReport(ME)
        fprintf(['******error in:runTest4mecka.m******\n']);
    end
end