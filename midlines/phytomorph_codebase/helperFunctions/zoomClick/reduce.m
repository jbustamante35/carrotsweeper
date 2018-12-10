function [dX,map] = reduce(X,b,map)

    if isempty(map)
        TL = [];
        ridx = randperm(size(X,4));
        for e = 1:10:min(size(X,4),8000)
            TL = [TL X(:,:,:,ridx(e))];
        end
        [iU,map] = rgb2ind(TL,2^b);
    end
    fprintf(['dithering data:']);
    for tr = 1:size(X,4)
        dX(:,:,tr) = dither(X(:,:,:,tr), map);
        if mod(tr,100) == 0
            fprintf(['.\n']);
        else
            fprintf(['.']);
        end
    end
    fprintf(['\n']);
    
end