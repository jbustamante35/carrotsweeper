function [fp] = filterP(image,psz)
    fp(1).image = image;
    for e = 1:numel(psz)
        tmpf = imfilter(image,
        domain = reshape(1:numel(image),size(image))
    end
end