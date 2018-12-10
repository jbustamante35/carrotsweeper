function [I] = func_depthStack(I,toT)
    % this function will permute (transpose) the image then
    % it will stack the color-panes [R G B] so that the 
    % column space spans the columns of the image
    if toT
        I = permute(I,[2 1 3]);
    end
    % permute before reshape
    I = permute(I,[1 3 2]);
    % get the size before reshape
    sz = size(I);
    % reshape the image
    I = reshape(I,[prod(sz(1:2)) prod(sz(3))]);
end