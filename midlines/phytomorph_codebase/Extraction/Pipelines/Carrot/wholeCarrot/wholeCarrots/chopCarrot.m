function [root top] = chopCarrot(MASK,chop,CONNECTLENGTH)
    % get the carrot top and root
    root = MASK(chop:end,:);    
    top = MASK(1:chop-CONNECTLENGTH,:);
    
    
    % check to see if there is a carrot
    R = regionprops(any(root,2),'Area','PixelIdxList');
    % find the index of the profile of the root
    fidx = find(sum(root,2));
end