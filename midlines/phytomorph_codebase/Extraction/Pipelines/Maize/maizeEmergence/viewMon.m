function [M] = viewMon(m,handF,algoF)
    if ischar(m)
        tmp = load(m);
        m = tmp.miniStack/255;
    end
    OFFSET = 4;
    H = cat(3,rgb2gray(m(:,:,:,handF-OFFSET)),rgb2gray(m(:,:,:,handF)),rgb2gray(m(:,:,:,handF+OFFSET)));
    A = cat(3,rgb2gray(m(:,:,:,algoF-OFFSET)),rgb2gray(m(:,:,:,algoF)),rgb2gray(m(:,:,:,algoF+OFFSET)));
    M = cat(1,H,A);
end