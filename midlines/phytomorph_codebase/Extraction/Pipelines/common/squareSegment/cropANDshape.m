function [d] = cropANDshape(inSpray,e0,e1,newYsize)
    try
        L = logical(any(inSpray{2},1));
        R = regionprops(L,'PixelIdxList');
        d = inSpray{1}(:,R(e1).PixelIdxList,:);
        d = imresize(d,[size(d,1) newYsize]);
        d = permute(d,[2 1 3]);
        % sort d
        d = sort(d,1);
        if ndims(d) == 3
            d = permute(d,[1 3 2]);
            sz = size(d);
            d = reshape(d,[sz(1)*sz(2) sz(3)]);
        end
    catch
        here = 1;
    end
end