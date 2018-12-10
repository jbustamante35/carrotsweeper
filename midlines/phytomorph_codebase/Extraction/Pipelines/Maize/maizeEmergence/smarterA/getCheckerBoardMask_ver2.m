function [MASK,boundingBox] = getCheckerBoardMask_ver2(I)
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the mask for the checker board
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['starting: getting checker board mask version 2.0\n']);tic;
        
        
        % hard filter on the red around the checkerboard
        targetColor = [220 120 80];
        targetHue1 = .04;
        HueWidth1 = .04;
        targetHue2 = .9;
        HueWidth2 = .2;
        
        
        
        % convert image to HSV
        hI = rgb2hsv(I);
        
        
        % filter based on target color values
        MASK = (abs(hI(:,:,1) - targetHue1) < HueWidth1 | ...
               abs(hI(:,:,1) - targetHue2) < HueWidth2) & hI(:,:,2) > .4;
        MASK = bwareaopen(MASK,100);
        MASK = bwlarge(MASK);
        MASK = imclose(MASK,strel('disk',11,0));
        MASK = imerode(MASK,strel('disk',11,0));
        R = regionprops(logical(MASK),'BoundingBox');
        boundingBox = R(1).BoundingBox;
        
        
        fprintf(['ending: getting checker board mask version 2.0: @ ' num2str(toc) '\n']);tic;
    catch ME
        fprintf(['failed on checkerboard mask \n']);
        getReport(ME)
    end
end