function [mask] = gatherPlantnHoods(I,mask,nhoodSZ,E,U,clusterG,n,clusterF)
    kI = [];
    fidx = find(mask);
    [ri ci] = find(mask);
    toSample = padarray(I,nhoodSZ,'replicate','both');
    KK2 = zeros(numel(ri),3*prod(2*(nhoodSZ)+1));
    %KK2 = zeros(numel(ri),6);
    cnt = 1;
    for p = 1:numel(ri)
        if cnt < size(KK2,1)
            stVec = toSample(ri(p)+nhoodSZ(1),(ci(p)+nhoodSZ(2)-nhoodSZ(2)):(ci(p)+nhoodSZ(2)+nhoodSZ(2)),:);
            stVec = sort(stVec,2);

            %{
            tmp = stVec;
            tmp = reshape(tmp,[1 2*nhoodSZ(2)+1 3]);
            f1 = mean(tmp(1,(end-4):end,:)) - mean(tmp(1,1:5,:));
            f2 = mean(tmp,2);
            f3 = [f1(:) f2(:)];
            %}
            KK2(cnt,:) = stVec(:);
            %KK2(cnt,:) = f3(:);
            cnt = cnt + 1;
        end
    end
    
    C = PCA_REPROJ(KK2,E,U);
    %{
    C = KK2;
 %}





    kidx = cluster(clusterG,C(:,n));
    
%{

    kidx2 = cluster(clusterF,C(kidx==3,4:6));

    kidx(kidx==3) = kidx2 + 2;

%}




    mask = double(mask);














    mask(fidx) = kidx;
end


