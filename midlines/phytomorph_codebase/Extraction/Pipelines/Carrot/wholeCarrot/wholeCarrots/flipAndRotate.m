function [I M] = flipAndRotate(I,M,EXP)
    if ~all(size(I) == EXP)
        I = permute(I,[2 1 3]);
        M = M';
    end
    tmp = imresize(M,.15);
    tmp = entropyfilt(tmp);
    tp = mean(mean(tmp(1:round(end/2),:)));
    bot = mean(mean(tmp(round(end/2):end,:)));
    if bot > tp
        M = flipud(M);
        I = flipdim(I,1);
    end
end