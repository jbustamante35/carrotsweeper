function [alignedBox] = alignBoxes(I,oPath,disp)
    
    % 1. Get the inside of the rectangle
    blackBoxMask = I(:,:,1) > 20 & I(:,:,1) < 105 & I(:,:,2) > 20 & I(:,:,2) < 105 & I(:,:,3) > 20 & I(:,:,3) < 105;
    softClose = imclose(blackBoxMask, strel('disk', 15, 0));
    filled = imfill(softClose, 'holes');
    I2 = filled == 1 & softClose == 0;
    I2 = imfill(I2, 'holes');

    % 2. Find each of the four corners
    [y,x] = find(I2);
    [~,loc] = min(y+x);
    C = [x(loc),y(loc)];
    [~,loc] = min(y-x);
    C(2,:) = [x(loc),y(loc)];
    [~,loc] = max(y+x);
    C(3,:) = [x(loc),y(loc)];
    [~,loc] = max(y-x);
    C(4,:) = [x(loc),y(loc)];

    % 3. Plot the corners
    % imshow(I); hold all
    % plot(C([1:4 1],1),C([1:4 1],2),'r','linewidth',3);

    %4. Find the locations of the new  corners
    L = mean(C([1 4],1));
    R = mean(C([2 3],1));
    U = mean(C([1 2],2));
    D = mean(C([3 4],2));
    C2 = [L U; R U; R D; L D];

    % 5. Do the image transform
    T = cp2tform(C,C2,'projective');
    IT = imtransform(I,T);
    I2R = imtransform(I2, T);
    out = flattenMaskOverlay(IT,logical(I2R)); %???
    I2R = bwlarge(I2R);

    % 6. Get rid of the lingering black
    I2R = imerode(I2R, strel('square', 21));
    R = regionprops(logical(I2R), 'boundingbox');
    alignedBox = imcrop(IT, R(1).BoundingBox);

end