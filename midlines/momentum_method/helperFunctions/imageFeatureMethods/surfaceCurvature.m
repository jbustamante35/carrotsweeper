function [centerP] = surfaceCurvature(imageStack,para)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute corner map
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %           I       := image
    %           para    := parameters for running the script         
    %                   := para.sig         -> sigma for gaussian filter
    %                   := para.gradPara    -> gradient configure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %           cM      := corner map       -> corner strength
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % I1 := image
    % I2 := MSK
    % I3 := objMsk    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % surface curvature
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute surface curvature
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear para;
    para{1} = 10;
    para{2} = .35;
    K = surKur(imageStack(:,:,1),para);
    deltaK = abs(K(:,:,1)-K(:,:,2));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute surface curvature local max points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear para;
    para.rad = 20;
    lMK = nonmaxsuppts(K(:,:,1), para.rad);
    lmK = nonmaxsuppts(K(:,:,2), para.rad);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % normalize curvature  maps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sK(:,:,1) = bindVec(K(:,:,1));
    sK(:,:,2) = bindVec(K(:,:,2));
    deltaK = bindVec(deltaK);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % logical call for surface curvature
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    para.minK_threshold = .8;
    logicMap = cat(3,sK(:,:,2) > para.minK_threshold,lmK);
    centerP = featureSelection(logicMap);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % surface curvature
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end