function [QR root] = splitBoxes(I,oPath,disp)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % start splitting the box
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tm = clock;
    fprintf(['Starting splitting the box' '\n']);


    redMask = I(:,:,1) > 115 & I(:,:,1) < 200 & I(:,:,2) > 70 & I(:,:,2) < 115 & I(:,:,3) > 95 & I(:,:,3) < 155;
    % Start filling in the borders (dilates and erode) 
    closeRedBox = imclose(redMask, strel('disk', 10, 0));
    % Fill in the boxes
    fillRedBox = imfill(closeRedBox, 'holes');
    % Get rid of the residuals pixels
    cleanRedBox = bwareaopen(fillRedBox, 5000);
    
    [rRed,cRed] = find(cleanRedBox);
    
    if disp
        imshow(cleanRedBox);
    end

    blueMask = I(:,:,1) > 53 & I(:,:,1) < 110 & I(:,:,2) > 130 & I(:,:,2) < 175 & I(:,:,3) > 170 & I(:,:,3) < 210;
    closeBlueCircle = imclose(blueMask, strel('disk', 11, 0));
    fillBlueCircle = imfill(closeBlueCircle, 'holes');
    cleanBlueCircle = bwareaopen(fillBlueCircle, 5000);
    
    [rBlue,cBlue] = find(cleanBlueCircle);
    
    if disp
        figure, imshow(cleanBlueCircle);
    end

    if (max(cRed) > min(cBlue))
        I = imrotate(I, 180);
        
        redMask = I(:,:,1) > 115 & I(:,:,1) < 200 & I(:,:,2) > 70 & I(:,:,2) < 115 & I(:,:,3) > 95 & I(:,:,3) < 155;
        closeRedBox = imclose(redMask, strel('disk', 10, 0));
        fillRedBox = imfill(closeRedBox, 'holes');
        cleanRedBox = bwareaopen(fillRedBox, 5000);
        [rRed,cRed] = find(cleanRedBox);
        
        blueMask = I(:,:,1) > 53 & I(:,:,1) < 110 & I(:,:,2) > 130 & I(:,:,2) < 175 & I(:,:,3) > 170 & I(:,:,3) < 210;
        closeBlueCircle = imclose(blueMask, strel('disk', 11, 0));
        fillBlueCircle = imfill(closeBlueCircle, 'holes');
        cleanBlueCircle = bwareaopen(fillBlueCircle, 5000);
        [rBlue,cBlue] = find(cleanBlueCircle);

    end
    
    splitPoint = mean([max(cRed), min(cBlue)]);
    QR = I(:,1:splitPoint,:);
    root = I(:,splitPoint:end,:);
    
    fprintf(['Box split: ' num2str(etime(clock,tm)) '\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % end splitting box
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%{
I = imread('/Users/boat/Dropbox/Wisconsin/Thesis/Phenotyping/image_preprocessing/2017/matlab_based/aligned/1_rep1_11.27.57_1.tif');
oPath = '/Users/boat/Desktop/output/';
[QR root] = splitBoxes(I, oPath, true);
%}
         
