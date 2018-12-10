function [I] = func_thumbNail(I,sz_per,toVec,toT)
    % function to resize the image and make thumbnail
    % if toT flag is true, it will transpose the data
    % if toVec flag is true, it will vectorize the data
    % resize the image
    I = imresize(I,sz_per,'nearest');
    % permute if needed
    if toT
        I = permute(I,[2 1 3]);
    end
    % vectorize if needed
    if toVec
        I = I(:);
    end
end