function [new_image] = mag(image,m)
    new_image = imresize(image,size(image)*m,'nearest');
end