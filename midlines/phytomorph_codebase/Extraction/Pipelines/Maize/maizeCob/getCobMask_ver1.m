function [Mask_out] = getCobMask_ver1(I,defaultAreaPix,colRange1,colRange2,fill)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                getCobMask_ver1.m turns image into black and white base upon background color information, 
                removes small objects(dusts) and operates closing.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                rgb2hsv_fast.m, imclearborder.m, imfill.m, bwareaopen.m, imopen.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                I:       An image to be analyzed in a matrix.
                defaultAreaPix: The default pixel to be considered noise relative to 1200 dpi.
                colRange1:      The color range for back ground to be removed in getcobMask.
                colRange2:      The color range for back ground to be removed in getcobMask.
                fill:           The radius of disk for Kernel of an image close operation.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    try
        fprintf(['start cob mask \n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % get gray and filter
        fprintf(['start rgb --> hsv conversion \n']);
        Mask_out = rgb2hsv_fast(I,'','H');
        fprintf(['end rgb --> hsv conversion \n']);
        Mask_out = Mask_out < colRange1/360 | Mask_out > colRange2/360;
        Mask_out = imfill(Mask_out,'holes');
        
        % changed for 300-1200 D.P.I on December 21, 2017
        % this might break the system, but rather than area-open
        % i will extract the X largest object
        %Mask_out = bwareaopen(Mask_out,defaultAreaPix);
        Mask_out = bwlarge(Mask_out,3);
        
        
        Mask_out = imopen(Mask_out,strel('disk',fill));
        % also removed the below line on December 21, 2017
        %Mask_out = bwareaopen(Mask_out,defaultAreaPix);
        fprintf(['end cob mask \n']);
    catch ME
        close all;
        getReport(ME);
        fprintf(['******error in:getCobMask_ver1.m******\n']);
    end
end