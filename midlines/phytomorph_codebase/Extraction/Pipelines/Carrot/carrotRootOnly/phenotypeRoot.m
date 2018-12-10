function [root] = phenotypeRoot(I,oPath,disp)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % start extracting profile
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tm = clock;
    fprintf(['Starting extracting the profile' '\n']);

    % Mask the blue circle
    blueMask = I(:,:,1) > 53 & I(:,:,1) < 110 & I(:,:,2) > 123 & I(:,:,2) < 175 & I(:,:,3) > 163 & I(:,:,3) < 210;
    closeBlueCircle = imclose(blueMask, strel('disk', 11, 0));
    fillBlueCircle = imfill(closeBlueCircle, 'holes');
    cleanBlueCircle = bwareaopen(fillBlueCircle, 10000);

    % Mask the background
    whiteMask = I(:,:,1) > 125 & I(:,:,2) > 125 & I(:,:,3) > 125;
    
    % Mask the object
    objectMask = cleanBlueCircle == 0 & whiteMask == 0;
    root = imopen(objectMask, strel('disk',5,0));
    root = imclearborder(root);
    root = bwareaopen(root, 500);
    
    % distortion = objectMask - cleanObject;
    % imshow(distortion);
    
    fprintf(['Profile extracted: ' num2str(etime(clock,tm)) '\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % end extracting profile 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

