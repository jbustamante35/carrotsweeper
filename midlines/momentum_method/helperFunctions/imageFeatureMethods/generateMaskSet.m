function [retM] = generateMaskSet(I)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make masks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make boundary mask
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear para;
    para.BW = [60 60 200 200];
    MSK = makeBoundaryMask(size(I),para.BW);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % subtract background from image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear para;
    para.closeVALUE = 31;
    para.resizeVALUE = .25;            
    bkI = subBackGrd(I,para);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % threshold, dilate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear para
    bkI = bindVec(bkI);
    level = graythresh(bkI);
    para.dilate = 31;
    objMsk = imdilate(bkI < level,strel('disk',para.dilate));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make masks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    retM.objMsk = objMsk;
    retM.bkI    = v;
    retM.MSK    = MSK;
end