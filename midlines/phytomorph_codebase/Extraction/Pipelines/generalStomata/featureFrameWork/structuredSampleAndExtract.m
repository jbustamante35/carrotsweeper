function [data] = structuredSampleAndExtract(rawI,toOp,trimValues,pointSet,opDomain,domainSize,extractionType,extractionArgs,disp)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rawI := raw image
    % toOp := image to operate on
    % trimValues := amount to snip off borders
    
    
    disp = false;
    if disp
        figH = figure;
    else
        figH = '';
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the inLine function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set download options
    options = weboptions('Timeout',60);
    % set local web path for saving loader program
    webPath = './cachePath/';
    mkdir(webPath);
    % set loader program
    inlineFeatureExtractionProgram = [webPath 'generalizeFeatureExtractor.mat'];
    % save generalized loader from irods
    websave(inlineFeatureExtractionProgram,'https://de.cyverse.org/dl/d/5C41DB00-9619-4FB0-A18B-66757FC1503A/generalizeFeatureExtractor.mat',options);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load the loader
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loader generalized loader
    featureP = load(inlineFeatureExtractionProgram);
    extractor = featureP.obj.func;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for r = 1:4
        rawI(:,1:trimValues(2)) = [];
        rawI = imrotate(rawI,90);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % peform extraction and analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['**********************************\n']);
    fprintf(['Starting opFunc over pointSet:\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for each slice
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for slice = 1:size(toOp,3)
        
        %tmpOp = toOp(:,:,slice);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for each point
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['start:Whole extracting bug eye image(s).\n']);
    str = clock;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subI = bugEyeViaT(toOp,pointSet(1,:),opDomain,domainSize,false);
    subI = thawTensor(subI,2);
    [T,NM] = extractor(subI,extractionType,extractionArgs);
    T = zeros(size(T,1),size(pointSet,1));

    %%%%%%%%%%%%%%%%%%
    parfor p = 1:size(pointSet,1)
        fprintf(['start:Extracting bug eye image ' num2str(p) ':' num2str(size(pointSet,1)) '.\n']);tic
        subI = bugEyeViaT(toOp,pointSet(p,:),opDomain,domainSize,disp,figH);
        subI = thawTensor(subI,2);
        [T(:,p)] = extractor(subI,extractionType,extractionArgs);
        eTime = toc;
        fprintf(['end:Extracting bug eye image ' num2str(p) ':' num2str(size(pointSet,1))  '@' num2str(eTime) '.\n']);
        rTime = mean(eTime)*(size(pointSet,1)-p);
        fprintf(['esti: time remaining : ' num2str(rTime) '\n']);
    end
    fprintf(['end:Whole extracting bug eye image(s) ' num2str(etime(clock,str)) '\n']);



    %T = opFunc(subI);
    %T = newF;
    %T = glueFrozenTensors(subI,newF);
        
    data.(NM) = T;
    data.rawI = rawI;
        
        
    fprintf('\n');
    fprintf(['Ending FFT:\n']);
    fprintf(['**********************************\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
