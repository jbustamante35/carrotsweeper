function [MASK] = getCheckerBoardMask(I)
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the mask for the checker board
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['starting: getting checker board mask \n']);tic;
        % hard filter on the red around the checkerboard
        targetColor = [220 120 80];
        targetHue = .04;
        HueWidth = .03;
        % convert image to HSV
        hI = rgb2hsv(I);
        % filter based on target color values
        MASK = abs(hI(:,:,1) - targetHue) < HueWidth & hI(:,:,2) > .4;
        MASK = bwareaopen(MASK,100);
        MASK = imclose(MASK,strel('disk',11,0));
        MASK = imerode(MASK,strel('disk',11,0));
        fprintf(['ending: getting checker board mask: @ ' num2str(toc) '\n']);tic;
    catch ME
        fprintf(['failed on checkerboard mask \n']);
        getReport(ME)
    end
end