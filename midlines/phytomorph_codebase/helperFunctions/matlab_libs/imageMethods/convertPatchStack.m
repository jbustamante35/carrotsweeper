function [tmp] = convertPatchStack(d)	
    for i = 1:size(d,ndims(d))
        tmp(:,:,i) = rgb2gray(d(:,:,:,i));
    end    
end