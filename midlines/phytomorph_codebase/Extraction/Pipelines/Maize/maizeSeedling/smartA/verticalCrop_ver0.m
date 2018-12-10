function [sI,sM,realS,R] = verticalCrop_ver0(I,MASK,newY)
    % get horizontal signals
    sig = any(MASK,1);
    % bulk up horizontal signal
    sig = repmat(sig,[size(I,1) 1]);
    % get regiongprops
    R = regionprops(sig);
    % for each region
    for e = 1:numel(R)
        % crop the image
        tmp = imcrop(I,R(e).BoundingBox);
        % get the orginal size of the horizontal direction
        realS(e) = size(tmp,2);
        % resize to the "choice" size
        tmp = imresize(tmp,[size(I,1) newY],'nearest');
        % store for return
        sI(:,:,:,e) = tmp;
        % crop the mask
        tmp = imcrop(MASK,R(e).BoundingBox);
        % resize to "choice" size
        tmp = imresize(tmp,[size(I,1) newY]);
        % store for return
        sM(:,:,e) = tmp;
    end
    
end