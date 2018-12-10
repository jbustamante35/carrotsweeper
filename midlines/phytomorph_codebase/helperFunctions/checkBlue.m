function [I] = checkBlue(I, checkBlue_scaleFactor,addcut,baselineBlue)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                checkBlue finds blue header from the raw image and if there is, it removes blue 
                header and returns resulted image.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                imresize.m, rgb2hsv_fast.m, regionprops.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:mecka('c',I800,3,oPutC,1,1,800);
                I:       An image to be analyzed in a matrix.
                checkBlue_scaleFactor:  A desired percentage to resize the image in checkBlue.
                rawImage_scaleFactor:   A desired percentage to resize the image.
                addcut:         The boarder handle for checkBlue. This is an addition to blue top computed in checkBlue.
                baselineBlue:   The baseline threshold to remove blue header in checkBlue.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    try
        fprintf(['starting with check for blue header.\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % hue and saturation
        % make mask
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        sI = imresize(I,checkBlue_scaleFactor);
        h = rgb2hsv_fast(sI,'single','H');
        v = rgb2hsv_fast(sI,'single','V');
        % find blue header from image 
        dI = h > 130/255 & h < 180/255 & v > 50/255;
        dI = sum(dI,2) > size(sI,2)/4;
        % measure blue area
        R = regionprops(dI,'Area');
        A = max([R.Area]);
        % modified September, 12 2017 - start
        if isempty(A)
            return;
        end
        % modified September, 12 2017 - end
        % it seems 600 would work down to 300DPI
        % it removes blue header from image if there exists
        fidx = find(dI);
        fidx = max(fidx);
        fidx = fidx *checkBlue_scaleFactor^-1;
        fidx = fidx + addcut;
        if A > baselineBlue*checkBlue_scaleFactor && fidx < size(I,1)/2;
            I = I(fidx:end,:,:);
        end
        fprintf(['end with check for blue header.\n']);
    catch ME
        getReport(ME)
        fprintf(['******error in checkblue.m******\n']);
    end
end