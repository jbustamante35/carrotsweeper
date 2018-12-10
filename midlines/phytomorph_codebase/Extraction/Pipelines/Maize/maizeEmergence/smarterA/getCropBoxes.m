function [boundingBox,centerPoints,MASK,I,uCL,covCL] = getCropBoxes(FileList,tform,topN,numToAverage,BOXSZ,cameraShift,disp,GMM)
    %{
    if isdeployed()
        parpool(1);
    end
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % transforming top N images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['starting: image transform of:' num2str(numToAverage)  '\n']);
    for s = 2:(numToAverage+1)
        fprintf(['starting rectification phase of:' num2str(s) ':' num2str(numToAverage) '\n']);tic
        if ~isempty(cameraShift)
            tmp = getRectifiedImage(FileList{s},tform,cameraShift{s});
        else
            tmp = getRectifiedImage(FileList{s},tform);
        end
        fprintf(['Starting RGB --> HSV \n']);
        rI(:,:,s-1) = rgb2hsv_fast(tmp,'','H');
        rIs(:,:,s-1) = rgb2hsv_fast(tmp,'','S');
        fprintf(['Ending RGB --> HSV \n']);
        fprintf(['ending rectification phase of:' num2str(s) ':' num2str(numToAverage) ':' num2str(toc) '\n'])
    end
    rI = mean(rI,3);
    rIs = mean(rIs,3);
    if ~isempty(cameraShift)
        I = getRectifiedImage(FileList{2},tform,cameraShift{1});
    else
        I = getRectifiedImage(FileList{2},tform);
    end
    
    fprintf(['ending: image transform of:' num2str(numToAverage)  '\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % transforming top N images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % hough transform
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['starting: hough transform \n']);
    OZ = round([95 130]);
    %{
    M = rI < .06;
    M = rI < .15 | rI > .7;
    %}
    
    % added on Sep 12 2018
    sz = size(I);
    M = reshape(I,[prod(sz(1:2)) sz(3)]);
    M = GMM.cluster(M);
    M = reshape(M,sz(1:2));
    M = M == 4;
   
    % removed on Sep 12 2018
    %M = rIs > .39 & (rI < .07 | rI > .16);
    M = bwareaopen(M,200);
    copyVec = M(end,:);
    M(end,:) = 0;
    M = imclearborder(M);
    M(end,:) = copyVec;
    M = imfill(M,'holes');
    %imwrite(M,'./output/whatdaf.tif');
    [centers, radii, metric] = imfindcircles(M,OZ,'ObjectPolarity','bright','Sensitivity',.99);
    mini = 4;
    CM = [];
    RM = [];
    for e = 1:size(centers,1)
        try
            CM(round(centers(e,2)/4),round(centers(e,1)/4)) = metric(e) + rand(1)*max(metric)*.0001;
            RM(round(centers(e,2)/4),round(centers(e,1)/4)) = radii(e);
        catch
            
        end
    end
    fprintf(['ending: hough transform \n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % hough transform
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % non max suppress circles based on hough strength
    % select out the top 168 circles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['starting: non-max suppression \n']);
    miniCM = imdilate(CM,strel('disk',round(1.75*mean(radii)/mini),0));
    %imwrite(miniCM,'./output/whatdaf2.tif');
    cp = (miniCM == CM).*(CM ~= 0);
    cp_index = find(cp);
    cp_strength = CM(cp_index);
    [J sidx] = sort(cp_strength,'descend');
    [cp1 cp2] = find(cp);
    %{
    if numel(cp1) < 100
        topN = topN / 2;
    end
    %}
    fprintf(['Reporting number of potential objects as:' num2str(numel(sidx)) '.\n']);
    cp_index = cp_index(sidx(1:topN));
    cp1 = cp1(sidx(1:topN));
    cp2 = cp2(sidx(1:topN));
    fprintf(['ending: non-max suppression \n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % non max suppress circles based on hough strength
    % select out the top 185 circles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make boundingBoxes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MASK = zeros(size(rI));
    for e = 1:size(cp1,1)
        boundingBox{e} = [mini*cp2(e) - BOXSZ  mini*cp1(e) - BOXSZ 2*BOXSZ 2*BOXSZ];
        MASK(round(mini*cp1(e)),round(mini*cp2(e))) = 1;
    end
    sampleMask = bwdist(MASK) < mean(.8*RM(cp_index));
    fidx = find(sampleMask);
    for k = 1:size(I,3)
        tmpI = I(:,:,k);
        CL(:,k) = tmpI(fidx);
    end
    uCL = mean(CL);
    covCL = cov(CL);
    MASK = imdilate(MASK,strel('disk',11,0));
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make boundingBoxes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % display
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if disp
        imshow(I,[]);
        viscircles(mini*[cp2 cp1],  RM(cp_index),'EdgeColor','r');
        for e = 1:numel(boundingBox)
            rectangle('Position',boundingBox{e},'EdgeColor','y')
        end
        drawnow
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % display
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    centerPoints = [cp2 cp1]*mini;
    
    rad_vec = 1.2*RM(cp_index);
end