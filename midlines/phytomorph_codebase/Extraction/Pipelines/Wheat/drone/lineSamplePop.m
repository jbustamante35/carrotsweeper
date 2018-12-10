function [F] = lineSamplePop(stack,F,sheets)
    tmp = zeros(size(sheets,1),size(sheets,2),size(stack,3));
    for k = 1:size(stack,3)
        tmp(:,:,k) = interp2(stack(:,:,k,end),sheets(:,:,1),sheets(:,:,2));
    end
    F(:,:,:,1) = [];
    
    F = cat(4,F,tmp);
end