function [vec tot] = frameMetric_ver0(fileName,gmm,N)
    
    tmp = imread(fileName);
    tot = rgb2gray(tmp);
    tot = sum(tot(:));
    tmp = rgb2hsv_fast(tmp,'','H');
    idx = cluster(gmm,tmp(:));
    idx = reshape(idx,[size(tmp,1) size(tmp,2)]);
    M = idx == N;
    M = bwareaopen(M,3000);
    M = imclearborder(M);
    vec = sum(M,2);
end