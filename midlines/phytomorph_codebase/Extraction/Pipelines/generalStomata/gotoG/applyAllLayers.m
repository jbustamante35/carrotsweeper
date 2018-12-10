function [data] = applyAllLayers(fileName,loaderType,loaderArgs,extractionType,extractionArgs,domainData,pointSetData,customData,AI_layer,basisU,basisE,optiPara,oPath,rPath)
    

    [pth,nm,ext] = fileparts(fileName);
    customData{1} = nm;


    [data] = extractionLayer(fileName,loaderType,loaderArgs,extractionType,extractionArgs,domainData,pointSetData,customData,[],[]);
    
    if isfield(data,'P')
        data.P = reshape(data.P',[432 432 size(data.P,1)]);
    end
    
    [~,BO,pMAP] = integrateAndmeasureProbMaps(data.P,[],optiPara,data.rawI);
    %{
    DT = bwdist(BO);
    sIDX1 = find(DT > 5);
    sIDX1 = randperm(numel(sIDX1));
    
    sIDX2 = find(DT > 50);
    sIDX2 = randperm(numel(sIDX2));
    
    %}
    
    
    cropImage = {};
    R = {};
    
    
    
    out = flattenMaskOverlay(data.rawI,logical(BO));
    
    EXPAND = 50;
    R = regionprops('table',BO,'BoundingBox','Centroid','Area');
    for r = size(R,1)
        R{r,'BoundingBox'}(1:2) = R{r,'BoundingBox'}(1:2) - EXPAND;
        R{r,'BoundingBox'}(3:4) = R{r,'BoundingBox'}(3:4) + 2*EXPAND;
        cropImage{r} = imcrop(data.rawI,R{r,'BoundingBox'});
    end
    
    
    cropImage = {};
    if ~isempty(oPath)
        % make the output directory
        mkdir(oPath);
        % save the mat file with struct option
        fileList{1} = [oPath customData{1} '.mat'];
        save(fileList{1},'cropImage','R','out','pMAP');
            
        % assign for push overlay
        fileList{end+1} = [oPath customData{1} '_overlay.jpg'];
        
        
        % write local
        imwrite(out,fileList{end});
        
        
        fileList{end+1} = [oPath customData{1} '_spotdata.csv'];
        writetable(R,fileList{end});
        
        
        
        % call push to irods
        pushToiRods(rPath,fileList);
    end
    
    
    %{
    tmpM = zeros(size(data.rI));
    tmpM(customData{2}) = 1;
    R = regionprops(logical(tmpM),'boundingBox');
    R(1).BoundingBox(3:4) = floor(R(1).BoundingBox(3:4)-1);
    subI = imcrop(data.rI,R(1).BoundingBox);
    %}
    % = applyAIlayer(data,AI_layer,basisU,basisE,[432 432]);
    
    
    %{
    init = [.1 1*3/9 1*3/9 1*3/9 200 .3 100 3 3];
    delta = [.05 .2 .2 .2 10 0.15 10 1 1]*.05^-1;
    mm{3} = @(X,Y)matthews_correlation(X,Y);
    mm{2} = @(X,Y)myMetric(X,Y);
    mm{1} = @(X,Y)positiveLikehoodRatio(X,Y);
    np = 5;
    paraLabels = [2*ones(1,np) 3*ones(1,np)];
    

    load('./fasterNET.mat','PT','para','PT2');

    for o = 1:numel(optiPara)
        [~,~,data.BO{o},data.overlayImage{o}] = integrateAndmeasureProbMaps(data.P,[],optiPara{o}(1:np),data.rawI,optiPara{o}((np+1):end),paraLabels,PT{o},PT2{o});
    end
    %}
    
    
    
    %{
    data.T = T;
    data.rI = rI;
    %}
    
    
end