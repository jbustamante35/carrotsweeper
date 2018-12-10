function [indexPosition] = generateImageDomain(tmpI,snip,skip)
    if nargin == 2
        skip = 1;
    end
    [d1,d2] = ndgrid((snip+1):skip:(size(tmpI,1)-snip),(snip+1):skip:(size(tmpI,2)-snip));
    % generate the index position(s) of the apply grid
    indexPosition = sub2ind(size(tmpI),d1(:),d2(:));
end