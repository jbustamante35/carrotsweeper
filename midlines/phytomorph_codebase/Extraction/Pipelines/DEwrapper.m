function [] = DEwrapper(functionName,varargin)
    fprintf(['First line of code in DEwrapper.\n']);
    [~,workingPath] = system('pwd');
    fprintf(['Working path is:' workingPath(1:end-1) '.\n']);
    date
    fprintf(['Version is: Debug-1.\n']);
    if toR('01-Apr-2019')
        OSG = false;
        selectedInput = varargin;
        
        fprintf(['*****************************\n']);
        fprintf(['Number of inputs:\n']);
        numel(selectedInput)
        fprintf(['*****************************\n']);
        fprintf(['start render of inputs:\n']);
        for e = 1:numel(selectedInput)
            selectedInput{e}
        end
        fprintf(['end render of inputs:\n']);
        fprintf(['*****************************\n']);
        
        
        fprintf(['****************************************************************************************\n']);
        versionString = ['Starting DE wrapper. \nPublication Version 1.0 - Sunday, April 3, 2017. \n'];
        fprintf(versionString);
        fprintf(['****************************************************************************************\n']);
        fprintf(['input1:' functionName '\n'])
        fprintf(['input2 Class:' class(selectedInput) '\n'])
        %{
        if nargin == 2
            OSG = false;
        end
        %}
        
        
        if isdir(selectedInput{1})
            filefunction = '';
            fileName = '';
            folderfunction = functionName;
            folderName = selectedInput{1};
        else
            filefunction = functionName;
            if numel(selectedInput) == 1
                fileName = selectedInput{1};
            else
                fileName = selectedInput;
            end
            folderfunction = '';
            folderName = '';
        end
        
        
        if ~isempty(filefunction)
            masterfuncToCall = filefunction;
        else
            masterfuncToCall = folderfunction;
        end
        
        %configurePeak(masterfuncToCall);
        
        %{
        fprintf(['input1:' filefunction '\n'])
        fprintf(['input2:' folderfunction '\n'])
        
        fprintf(['input3:' fileName '\n'])
        fprintf(['input4:' folderName '\n'])
        %}
        
        if isempty(filefunction) && isempty(folderfunction)
            fprintf(['You must select only one algorithm at a time \n']);
            return
        end
        
        
        if ~isempty(filefunction)
            fprintf(['running:filefunction:' filefunction '@' num2str(numel(fileName)) '\n']);
        end
        
        
        if ~isempty(folderfunction)
            fprintf(['running:folderfunction:' folderfunction '@' folderName '\n']);
        end
    
        
        pwd = getenv('PWD');
        fprintf(['Current working dir:' pwd '\n']);
        
        options = weboptions('Timeout',60);

        
        
        
        webPath = './';
        if ~OSG
            webPath = '/de-app-work/';
            sPATH = '//usr/bin/';
            
            
            % set for dcraw
            fprintf(['start:setting up dcraw:\n']);
            fprintf(['start:setting up dcraw:1:downloading dcraw\n']);
            websave([sPATH 'dcraw'],'http://de.cyverse.org/dl/d/B6A77367-7921-42EC-ADAC-1C44EEE51331/dcraw',options);
            CMD = 'chmod +x /usr/bin/dcraw';
            [r,q] = system(CMD,'-echo');

            % setup for libs for dcraw
            fprintf(['start:setting up dcraw:2:setting up libs\n']);
            websave('/usr/bin/lcms_lib.tar.gz','http://de.cyverse.org/dl/d/24C561A0-D452-4A23-BCE2-311468A4867F/lcms_lib.tar.gz',options);
            CMD = 'tar xvf /usr/bin/lcms_lib.tar.gz -C /usr/bin/';
            [r,q] = system(CMD,'-echo');

            % setup ld_library_path
            fprintf(['start:setting up dcraw:3:setting up LD_LIBRARY_PATH\n']);
            OLD_PATH = getenv('LD_LIBRARY_PATH');
            setenv('LD_LIBRARY_PATH', ['/usr/lib/x86_64-linux-gnu/:' getenv('LD_LIBRARY_PATH') ':/usr/bin/lcms/lib/']);
            [r,q] = system(CMD,'-echo');
            fprintf(['LD_LIBRARY_PATH IS:\n']);
            CMD = 'echo $LD_LIBRARY_PATH';
            [r,q] = system(CMD,'-echo');
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % setup libjpeg62
            CMD = 'apt-get install libjpeg62';
            [r,q] = system(CMD,'-echo');

            % setup java libs for QR
            fprintf(['start:setting up javase-3.2.1.jar:1:downloading javase-3.2.1.jar\n']);
            websave('/de-app-work/javase-3.2.1.jar','http://de.cyverse.org/dl/d/385C747A-C33D-40AB-BEF1-B62279FA8785/javase-3.2.1.jar',options);

            % setup java libs for QR
            fprintf(['start:setting up core-3.2.1.jar:1:downloading core-3.2.1.jar\n']);
            websave('/de-app-work/core-3.2.1.jar','http://de.cyverse.org/dl/d/10CC5811-B22A-44CF-87DC-BC4464A8D5BB/core-3.2.1.jar',options);

            % setup java libs for QR
            fprintf(['start:setting up core-3.2.1.jar:1:downloading bioformats_package.jar\n']);
            websave('/de-app-work/bioformats_package.jar','http://de.cyverse.org/dl/d/97063E78-CF12-489F-9BFB-1FE01C71912A/bioformats_package.jar',options);
            
            % setup Cory Hirsch perl script for QR gen
            fprintf(['start:perl script for QR:1:generate_qrcode_R_script.pl\n']);
            websave('/de-app-work/generate_qrcode_R_script.pl','https://de.cyverse.org/dl/d/3DB5C874-623A-41F6-8108-05D183C77B42/generate_qrcode_R_script.pl',options);
            CMD = ['chmod +x /de-app-work/generate_qrcode_R_script.pl'];
            [status,cmdout] = system(CMD,'-echo');
            
            % restore path
            setenv('LD_LIBRARY_PATH', [OLD_PATH ':/usr/bin/lcms/lib/']);
            fprintf(['LD_LIBRARY_PATH IS:\n']);
            CMD = 'echo $LD_LIBRARY_PATH';
            [r,q] = system(CMD,'-echo');
            
            
            if isdeployed
                javaaddpath([pwd filesep 'core-3.2.1.jar']);
                javaaddpath([pwd filesep 'javase-3.2.1.jar']);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
        
        
        
        
        
        
        if ~isempty(folderfunction)
            % dig for images
            folderName = [folderName filesep];
            FileList = {};
            FileExt = {'tif','TIF','BMP','bmp','PNG','png','JPG','jpg','tiff','csv','json'};
            FileList = gdig(folderName,FileList,FileExt,1);
            toOpOn = FileList;
            switch folderfunction
                case 'midResRoot'
                    % init values
                    outputLocation = './output/';
                    toWrite = '1';
                    NP = '20';
                    NPK = '20';
                    SNIP = '20';
                    disp = '1';
                    [pth,nm,ext] = fileparts(FileList{1});
                    pth = strrep(pth,filesep,'SLASH');
                    matName = [pth nm];
                    % create anon funcToCalltion
                    funcToCall = @(X)isolateRoots_overStack(cellStringConvert(X),outputLocation,matName,toWrite,NP,NPK,SNIP,disp,'-1');
                case 'hypocotylSegments'
                    outputLocation = './output/';
                    disp = 0;
                    funcToCall = @(X)op0_liteDE(X,outputLocation,disp);
                case 'rootPH'
                    outputLocation = './output/';
                    funcToCall = @(X)singlePHstack(X,outputLocation);
                case 'wholeCarrot-Stage2'
                    fprintf(['start:pulling petLength.mat\n']);
                    websave([webPath 'petLength.mat'],'https://de.cyverse.org/dl/d/615CA1A4-CD2F-4F25-9E8E-0AA0BDCD9639/petLength.mat',options);
                    load([webPath 'petLength.mat']);
                    funcToCall = obj.func;
                    delete([webPath 'petLength.mat']);
                case 'TIPS'
                    funcToCall = @(X)folderTIPS(X);
                case 'imageClassifier'
                    funcToCall = @(X)generateImageClassifier(X);
                case 'maizeGerm'
                    fprintf(['start:pulling cornPopperNNapp.mat\n']);
                    %websave([webPath 'cornPopperNNapp.mat'],'https://de.cyverse.org/dl/d/E396B8C2-C059-4405-879B-D8B1095DD89A/cornPopperNNapp.mat',options);
                    websave([webPath 'cornPopperNNapp.mat'],'https://de.cyverse.org/dl/d/6EB4ADC9-1D8F-424E-95D5-03E2C715B142/cornPopperNNapp.mat',options);
                    load([webPath 'cornPopperNNapp.mat']);
                    funcToCall = obj.func;
                    delete([webPath 'cornPopperNNapp.mat']);
                case 'maizeSwell'
                    funcToCall = @(X)measureStack_finalDeploy(X,[200 200 200 200],'./output/',false);
                case 'overHeadCamera'
                    %{
                    fprintf(['start:pulling overHeadfuncToCallAPP.mat\n']);
                    websave([webPath 'overHeadfuncToCallAPP.mat'],'https://de.cyverse.org/dl/d/77995EAD-52E8-446D-867A-631D2BCA5339/overHeadfuncToCallAPP.mat',options);
                    load([webPath 'overHeadfuncToCallAPP.mat']);
                    funcToCall = obj.func;
                    delete([webPath 'overHeadfuncToCallAPP.mat']);
                    %}
                    fprintf(['start:pulling overHeadfuncToCallAPP.mat\n']);
                    websave([webPath 'overHeadfuncToCallAPP.mat'],'https://de.cyverse.org/dl/d/7507C2DA-9F2F-41AB-8408-888670320E71/overHeadFuncAPP2.mat',options);
                    load([webPath 'overHeadfuncToCallAPP.mat']);
                    funcToCall = obj.func;
                    delete([webPath 'overHeadfuncToCallAPP.mat']);
                case 'kernelScanalyzer'
                    fprintf(['start:pulling scanalyzer.mat\n']);
                    websave([webPath 'scanalyzer.mat'],'https://de.cyverse.org/dl/d/3678726E-63E4-4AE2-BA93-D9D8B012E055/scanalyzer.mat',options);
                    load([webPath 'scanalyzer.mat']);
                    funcToCall = obj.func;
                    delete([webPath 'scanalyzer.mat']);
            end
        end
        
        
        
        if ~isempty(filefunction) 
            toOpOn = fileName;
            switch filefunction
                case 'featurePoints'
                    funcToCall = @(X)masterSpore(X,60,'./output/');
                case 'arabidopsisForDirk'
                    funcToCall = @(X)cellForDirk(X,'./output/');
                case 'JSONcompile'
                    funcToCall = @(X)JSONcompile(X,'./output');
                case 'QRTEXTcompile'
                    funcToCall = @(X)QRcompile(X,'./output');
                case 'sorghumLeaf'
                    funcToCall = @(X)sorguhmLeafAnalysis(X,'./output/');
                case 'confocalRootTip'
                    funcToCall = @(X)confocalRatio(X,1,'./output/',0);
                case 'wholeCarrot'
                    funcToCall = @(X)singleWholeCarrotAnalyze(X,3350,200,40,'./output/',[]);
                case 'singleEarImage'
                    toOpOn = fileName;
                    algorithm = 'e';
                    numberOfObjects = '3';
                    oPath = './output/';
                    remotePath = '';
                    toSave = '1';
                    toDisplay = '1';
                    scanResolution = '1200';
                    rawImage_scaleFactor = 1;
                    funcToCall = @(X)mecka(algorithm,X,numberOfObjects,oPath,remotePath,toSave,toDisplay,scanResolution,rawImage_scaleFactor);
                case 'singleCobImage'
                    algorithm = 'c';
                    numberOfObjects = '3';
                    oPath = './output/';
                    remotePath = '';
                    toSave = '1';
                    toDisplay = '1';
                    scanResolution = '1200';
                    rawImage_scaleFactor = 1;
                    funcToCall = @(X)mecka(algorithm,X,numberOfObjects,oPath,remotePath,toSave,toDisplay,scanResolution,rawImage_scaleFactor);
                case 'singleKernelImage'
                    algorithm = 'k';
                    numberOfObjects = [];
                    oPath = './output/';
                    remotePath = '';
                    toSave = '1';
                    toDisplay = '1';
                    scanResolution = '1200';
                    rawImage_scaleFactor = 1;
                    funcToCall = @(X)mecka(algorithm,X,numberOfObjects,oPath,remotePath,toSave,toDisplay,scanResolution,rawImage_scaleFactor);
                case 'maizeSeedling'
                    %{
                    smoothValue = 100;
                    threshSIG = 5;
                    EXT = 400;
                    topTRIM = 100;
                    SNIP = 4;
                    BKBOUND = 20;
                    outputLocation = './output/';
                    remoteOutputLocation = [];
                    funcToCall = @(X)singleSeedlingImage(X,smoothValue,threshSIG,EXT,topTRIM,SNIP,BKBOUND,outputLocation,remoteOutputLocation);
                    %}
                    %{
                    fprintf(['start:pulling maizeSeedlings.mat\n']);
                    websave([webPath 'maizeSeedlings.mat'],'https://de.cyverse.org/dl/d/7B2CBBBA-EA62-4DA6-836E-8093366BD8BE/maizeSeedlings.mat',options);
                    load([webPath 'maizeSeedlings.mat']);
                    %}
                    
                    fprintf(['start:pulling maizeSeedlings_HOTFIX.mat\n']);
                    websave([webPath 'maizeSeedlings_HOTFIX.mat'],'https://de.cyverse.org/dl/d/69389EB9-FDC1-4A64-BE4B-D975E2CE0634/maizeSeedlings_HOTFIX.mat',options);
                    load([webPath 'maizeSeedlings_HOTFIX.mat']);
                    
                    funcToCall = obj.func;
                    delete([webPath 'maizeSeedlings_HOTFIX.mat']);
        
                case 'sorghumStomata'
                    
                    fprintf(['start:pulling sorghumNNapp.mat\n']);
                    websave([webPath 'sorghumNNapp.mat'],'https://de.cyverse.org/dl/d/A2E6EFA2-37C3-4EC6-9335-302C29F9B2CC/metaLayers_sorghum.mat',options);
                    load([webPath 'sorghumNNapp.mat']);
                    funcToCall = obj.func;
                    delete([webPath 'sorghumNNapp.mat']);
                    
                    %{
                    fprintf(['start:pulling sorghumNNapp.mat\n']);
                    websave([webPath 'sorghumNNapp.mat'],'http://de.cyverse.org/dl/d/D9A3E7EB-7260-4E1E-BE44-8DF8BD0FA08B/sorghumNNapp.mat',options);
                    load([webPath 'sorghumNNapp.mat']);
                    funcToCall = obj.func;
                    delete([webPath 'sorghumNNapp.mat']);
                    %}
                    
                    %{
                    [pth,nm,ext] = fileparts(toOpOn);
                    customData{1} = ['{name_' strrep(pth,filesep,'FILESEP') 'FILESEP' nm '}{patch_' num2str(1) '}'];
                    funcToCall = @(X)cRunner(X,'./','',customData);
                    %}
                case 'maizeStomata'
                    %{
                    fprintf(['start:pulling maizeStomataNNapp.mat\n']);
                    websave([webPath 'maizeStomataNNapp.mat'],'https://de.cyverse.org/dl/d/515B2440-3221-4D60-98D8-6FD4DCEA7CA7/maizeStomataNNapp.mat',options);
                    load([webPath 'maizeStomataNNapp.mat']);
                    funcToCall = obj.func;
                    delete([webPath 'maizeStomataNNapp.mat']);
                    %}
                    fprintf(['start:pulling maizeStomataNNapp.mat\n']);
                    websave([webPath 'maizeStomataNNapp.mat'],' https://de.cyverse.org/dl/d/4FCB6C9C-93E3-49C3-BE85-47AB1854C0CB/metaLayers_corn.mat',options);
                    load([webPath 'maizeStomataNNapp.mat']);
                    funcToCall = obj.func;
                    delete([webPath 'maizeStomataNNapp.mat']);
                   
                case 'arabidopsisSeed'
                    funcToCall = @(X)measureArabidopsisSeed(X,'./output/');
                case 'DNAcrossover'
                    funcToCall = @(X)measureCrossOver(X,'./output/','');
                case 'generateQRlabels'
                    funcToCall = @(X)genQRsheets_latex(X,'./output/');
                case 'generateQRsheets'
                    funcToCall = @(X)genQRLargeFormatSheets(X,'./output/');
            end 
        end
        
        
        
        if ~isempty(filefunction)
            fprintf(['running:filefunctiontion:' filefunction '@' num2str(numel(fileName)) '\n']);
        end
        
        
        if ~isempty(folderfunction)
            fprintf(['running:folderfunction:' folderfunction '@' folderName '\n']);
        end
    
        
        fprintf(['Start from DE.\n']);
        funcToCall(toOpOn);
        fprintf(['End from DE.\n']);
        fprintf(['Cleaning up files from disk.\n']);
        % clean up
        delete([webPath 'generate_qrcode_R_script.pl'],[webPath 'javase-3.2.1.jar'],[webPath 'core-3.2.1.jar'],[webPath 'bioformats_package.jar']);
        fprintf(['Last line of code.\n']);
    end
    
    
end

%{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TEST
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %-a SeriesNetwork.m
    funcToCalltionFile = '-a im2java2d.m -a detectEmergence_ver2.m -a rgbAndhsvExtract.m -a arborist.m -a measureCrossOver.m -a extractPhenotypesFromOverheadCamera_ver2.m -a clusterForest.m -a clusterTree.m -a clusterNode.m -a extractPhenotypesFromOverheadCamera.m -a genQRLargeFormatSheets.m -a detectEmergence -a imsubtract.m -a rgb2hsv_fast.m -a generateImageClass.m -a getPNN_funcToCall.m -a shapeVerticalStripNozzle.m -a network.m -a funcToCall_depthStack.m -a funcToCall_resizeDepthStack.m -a constantTransitionfuncToCalltion.m -a myProb.m -a my_hmm.m -a hmm_node.m -a nozzleManifold.m -a cropImages_v2.m -a funcToCall_thumbNail.m -a smartMain_v4.m -a smartMain_v2.m -a sorguhmLeafAnalysis.m -a maizeSeedling_funcToCall3.m -a maizeSeedling_funcToCall2.m -a maizeSeedling_funcToCall1.m -a smartMain.m -a singleWholeCarrotStage2.m -a petLength.m -a partialfuncToCalltion.m -a condorDeploy_ver0.m -a confocalRatio.m -a isolateRoots_overStack.m -a mecka.m -a singleWholeCarrotAnalyze.m -a op0.m -a singleSeedlingImage.m';
    cdir = dir('/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/Maize/maizeSeedling/nets/');
    cdir(1:2) = [];
    for e = 1:numel(cdir)
        funcToCalltionFile = ['-a ' cdir(e).name ' ' funcToCalltionFile];
    end
    tmpCompileDirectory = '/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/deDeployR2017a/';
    mkdir(tmpCompileDirectory)
    tmpCompileDirectory = '/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/deDeploy/';
    mkdir(tmpCompileDirectory)

    CMD = ['mcc -d ' tmpCompileDirectory ' ' funcToCalltionFile ' -a  gmdistribution.m -a cJob.m -a im2single.m -m -v DEwrapper.m'];
    eval(CMD);
    pushCMD = ['iput -f /mnt/scratch1/phytomorph_dev/Extraction/Pipelines/deDeploy/DEwrapper /iplant/home/nmiller/publicData/DEwrapper'];
    [pushR] = system(pushCMD,'-echo');
    iput from here:/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/deDeploy
    iput to here:/iplant/home/nmiller/publicData

    W = functions(obj.func);
    W.workspace


%}