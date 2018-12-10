function [] = runPipelines()
    PipelineList = {'Arabidopsis Mid Resolution','Maize','Maize Swelling'};   
    pipelineSelection = printList(PipelineList,'Pipelines','Please select a pipeline to run:');
    
    switch pipelineSelection
        case 1
            runArabidopsisPipeline();
        case 2
            runMaizePipeline();
        case 3
            runSwellingPipeline();
    end
end

function [selection] = printList(List,Title,Prompt)
    fprintf([Title '\n']);
    for e = 1:numel(List)
        fprintf(['[' num2str(e) ']-' List{e} '\n']);
    end
    selection = input(Prompt);
end

function [] = runArabidopsisPipeline()
    in.inScanFile = '/mnt/snapper/pipeLines/arabidopsis/scanPaths.csv';
    in.processedFile = '/mnt/snapper/pipeLines/arabidopsis/processed.csv';
    scanForImages(in);
end

function [] = runMaizePipeline()
    MaizePipelineList = {'Maize Ears','Maize Cobs','Maize Kernels','All','None'};
    MaizeUserList = {'nhaase','gxe','ruairidh','garf0012','All'};
    FunctionList = {'Sync Data','Run Algorithm','Sync and Run Algorithm'};
    user = printList(MaizeUserList,'User List','Please select users pipline:');
    method = printList(MaizePipelineList,'Method List','Please select algorithm:');
    func = printList(FunctionList,'Function List','Please select function:');
    if method == 4;method = 1:3;end
    if user == numel(MaizeUserList)+1;user= 1:numel(MaizeUserList);end
    if func == 3;func == 1:2;end
    
    for u = 1:numel(user)
        curUser = MaizeUserList{user(u)};
        for f = 1:numel(func)
            curFunc = FunctionList{func(f)};
            switch curFunc
                case FunctionList{1} % sync data
                    fprintf('No longer Supported');                    
                case FunctionList{2} % run algorithm
                    for m = 1:numel(method)
                        curMeth = method(m);
                        switch curMeth
                            case 1
                                scanAndAnalyzeMaizeEars(curUser)
                            case 2
                                scanAndAnalyzeMaizeCobs(curUser);
                            case 3
                                scanAndAnalyzeMaizeKernels(curUser);
                        end     
                    end
            end
            
        end
    end
end

function [] = runSwellingPipeline()
    SwellingProjectList = {'Scott','Het','IBM_syn10','Cali','Gusheng'};
    project = printList(SwellingProjectList,'Project List','Please project pipeline:');
    % raw path for input and output paths
    inFilePath = '/mnt/snapper/kernelSwellingData/#project#/rawData/';
    %oPath = '/mnt/snapper/kernelSwellingData/#project#/return_second_WIDIV/';
    oPath = '/mnt/snapper/kernelSwellingData/#project#/return/';
    % replace the #project#
    inFilePath = strrep(inFilePath,'#project#',SwellingProjectList{project});
    oPath = strrep(oPath,'#project#',SwellingProjectList{project});
    mkdir(oPath);
    
    % this is the run for the first dataset
    %main_swell_HT(inFilePath,oPath,9,95,95,@(x,y)lt(x,y));
    template1 = [];
    template2 = [];
    for r = 1:10
        template1 = [template1;r*ones(1,4)];
    end
    MIDDLE = 11*ones(10,1);
    for r = 12:21
        template2 = [template2;r*ones(1,4)];
    end
    template = [template1 MIDDLE template2];
    template1 = [];
    for r = 1:10
        template1 = [template1;r*ones(1,8)];
    end
    template = template1;
    template = [];
    for r = 1:10
        template = [template;r*ones(1,9)];
    end
    % this is the run for the big proper datasets
    %ExpectedNumberOfImages = 144;
    %toRun = 144;
    ExpectedNumberOfImages = 130;
    toRun = 130;
    kernelsPerRow = 9;
    main_swell_HT(inFilePath,oPath,kernelsPerRow,ExpectedNumberOfImages,toRun,@(x,y)lt(x,y),1,template);
end