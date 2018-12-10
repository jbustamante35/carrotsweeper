function [geoP] = geodesicCurvature(imageStack)
    
    % I1 := image
    % I2 := MSK
    % I3 := objMsk    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % geodesic curvature
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute surface curvature
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear para;
    para{1} = 5;
    para{2} = .85;
    gK = geoKur(imageStack(:,:,1),para);
    
    expK = mean(gK(:));
    gK = imgReplace(gK,~imageStack(:,:,2),expK);
    gK = imgReplace(gK,gK < 0,expK);
    gK = bindVec(gK);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute surface curvature local max points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear para;
    para.rad = 7;
    lMgK = nonmaxsuppts(gK, para.rad);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % logical call for geodesic curvature
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear para
    para.minK_threshold = .1;
    logicMap = cat(3,gK > para.minK_threshold,lMgK,imageStack(:,:,2:end));
    geoP = featureSelection(logicMap);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % geodesic curvature
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end