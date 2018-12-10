function [M spec] = makeKernelMask_ver0(image)
    ng = 3;
    sz = size(image);
    F = reshape(image,[prod(sz(1:2)) sz(3)]);
    [S C U E L ERR LAM] = PCA_FIT_FULL(F,3);    
    BK = kmeans(C,ng);
    
    BK = reshape(BK,sz(1:2));
    bkidx = mode(BK(end,:));
    %{
    stack = [];
    for l = 1:ng
        stack = [stack sum(BK==l)];
    end
    [J bkidx] = max(stack);
    BK = reshape(BK,sz(1:2));
    %}
    M = ~(BK == bkidx);
    %{
    
    mM = imerode(M,strel('disk',21,0));
    R = regionprops(M,'Centroid','PixelIdxList','Area');    
    cidx = count([R.Area]);
    
    R = R(cidx==1);
    %Z = zeros(size(BK));
    for r = 1:numel(R)
        %{
        Z = zeros(size(BK));
        Z(R(r).PixelIdxList) = 1;        
        imshow(Z,[])
        drawnow
        %}
        spec.dataS(r,:) = std(F(R(r).PixelIdxList,:),1,1);
        spec.data(r,:) = mean(F(R(r).PixelIdxList,:),1);
        spec.loc(r,:) = R(r).Centroid;
    end
    %}
    spec = [];
end