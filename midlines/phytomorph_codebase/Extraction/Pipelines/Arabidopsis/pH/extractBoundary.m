function [dB] = extractBoundary(fileName,filtNumber)
    I = imread(fileName);
    I = I(:,:,1:3);
    I = permute(I,[2 1 3]);
    G = rgb2gray(I);
    M = isolateRoot3((I(:,:,2)),@(x)entropyfilt(x,ones(filtNumber)),filtNumber,0);
    dB = bwboundaries(M);
    rmidx = find(dB{1}(:,2) == 1);
    dB{1}(rmidx,:) = [];
end