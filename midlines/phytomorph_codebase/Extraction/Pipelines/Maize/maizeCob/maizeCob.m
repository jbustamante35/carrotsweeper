function [BB PS] = maizeCob(I,numberOfCobs,defaultAreaPix,colRange1,colRange2,fill)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                maizeCob.m 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                getCobMask_ver1.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                I:       An image to be analyzed in a matrix.
                numberOfCobs:            Number of cobs that are expected to be analyzed. 
                defaultAreaPix: The default pixel to be considered noise relative to 1200 dpi.
                colRange1:      The color range for back ground to be removed in getcobMask.
                colRange2:      The color range for back ground to be removed in getcobMask.
                fill:           The radius of disk for Kernel of an image close operation.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    %%%%%%%%%%%%%%%%%%%%%%%
    % init return vars    
    BB = [];
    PS = [];
    %%%%%%%%%%%%%%%%%%%%%%%
    try
        fprintf(['starting single maizeCob.m analysis \n']);
        % declare vars
        PS.widthProfile = [];
        % get cob mask
        [B] = getCobMask_ver1(I,defaultAreaPix,colRange1,colRange2,fill);
        R = regionprops(B,'PixelIdxList','PixelList','Area','Image','BoundingBox');
        if numel(R) > numberOfCobs
            [~,sidx] = sort([R.Area]);
            R(sidx(1:(end-numberOfCobs))) = [];
        end
        numberOfCobs = min(numberOfCobs,numel(R));
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % iterate over the cobs    
        for e = 1:numberOfCobs
            fprintf(['starting cob analysis for object:' num2str(e) '\n']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % get bounding box and the image
            BB{e} = R(e).BoundingBox;
            % get the binary mask
            subI = R(e).Image;
            
            % crop the color image for obtaining the cob color
            tmpI = imcrop(I,floor((R(e).BoundingBox)));            
            % make strip image for cob color
            tmpM = zeros(size(subI));
            tH = R(e).BoundingBox(4)/3;            
            tStr = round(+tH);
            tStp = round(2*tH);
            tmpM(tStr:tStp,:) = 1;
            % fill in the middle third
            tmpM = tmpM.*subI;
            PixelIdxList = find(tmpM);
            for k = 1:size(I,3)
                tmpP = tmpI(:,:,k);
                tmpC(k) = mean(tmpP(PixelIdxList));
            end
            RGB(e,:) = tmpC;
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % get width profile
            WIDTH = sum(subI,2);
            fidx = find(WIDTH ~= 0);
            uW = mean(WIDTH(fidx));
            PS.average_WIDTH(e) = uW;
            PS.widthProfile = [PS.widthProfile ;interp1(1:numel(WIDTH),WIDTH,linspace(1,numel(WIDTH),1000))];
        end
        PS.RGB = RGB;
    catch ME
        close all;
        getReport(ME);
        fprintf(['******error in:maizeCob.m******\n']);
    end
end
