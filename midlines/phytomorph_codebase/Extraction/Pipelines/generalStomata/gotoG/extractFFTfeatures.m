function [T] = extractFFTfeatures(F,NF,N)
    F = thawTensor(F,2);
    M = mfftm(N(2),NF);
    F = permute(F,[2 1 3]);
    fsz = size(F);
    fF = mtimesx(M,F);
    fF = reshape(fF,[size(M,1) fsz(2:3)]);
    fF = permute(fF,[2 1 3]);
    T{1} = abs(fF);
    T{2} = angle(fF);
    T = freezeTensor(T,true);
end