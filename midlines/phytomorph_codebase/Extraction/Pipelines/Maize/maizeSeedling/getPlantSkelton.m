function [skeleton] = getPlantSkelton(MASK)
    % pad the array
    tmpMASK = padarray(MASK, [300 0], 'replicate', 'post');
    % get the skeleton
    skeleton = bwmorph(tmpMASK,'thin',inf);
    % get the skeleton
    skeleton = skeleton(1:size(MASK,1),:);
    skeleton = bwareaopen(skeleton,10);
    skeleton = bwlarge(skeleton);
end