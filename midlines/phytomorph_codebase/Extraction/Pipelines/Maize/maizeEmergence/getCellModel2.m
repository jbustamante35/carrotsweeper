function [idx,model,POS,out] = getCellModel2(I,obj,obj2,capIDX,spaceClusterFunc,disp)
    sz = size(I);
    zMASK = I(:,:,1) ~= 0;
    out = I;
    I = reshape(I,[prod(sz(1:2)) sz(3)]);
    [idx,~,prob] = cluster(obj,I);
    fidx1 = find(idx == capIDX);






    fidx2 = find(idx ~= capIDX);
    [idx2,~,prob2] = cluster(obj2,I(fidx2,:));

    idx(fidx1) = 1;
    idx(fidx2) = idx2 + 1;
    
    newPROB = zeros(size(prob,1),size(prob2,2));
    newPROB(fidx2,:) = prob2;
    PROB = [prob(:,capIDX) newPROB];
    prob = PROB;


    idx = reshape(idx,sz(1:2));
    idx = idx.*double(zMASK);
    CL = {'b' 'r' 'b'};
    
    for k = 1:3
        out(:,:,k) = bindVec(out(:,:,k));
    end
    
    [n1 n2] = ndgrid(1:size(idx,1),1:size(idx,2));
    
    
    for k = 1:size(prob,2)
        model(:,:,k) = reshape(prob(:,k),sz(1:2));
        model(:,:,k) = model(:,:,k).*zMASK;
        cl(:,:,k) = spaceClusterFunc{k}(model(:,:,k));
        cl(:,:,k) = cl(:,:,k) / sum(sum(cl(:,:,k)));
        tmpP = cl(:,:,k);
        POS(1,k) = n1(:)'*tmpP(:);
        POS(2,k) = n2(:)'*tmpP(:);
        if disp
            msk = zMASK & idx == k;
            %msk = bwlarge(msk);
            out = flattenMaskOverlay(out,msk,.3,CL{k});
        end
    end

    

    

%{
    model(:,:,1) = reshape(mvnpdf(I,uCap,sCap),sz(1:2));
    model(:,:,2) = reshape(mvnpdf(I,uPVC,sPVC),sz(1:2));
%}
    %model(:,:,1) = bindVec(model(:,:,1));
    %model(:,:,2) = bindVec(model(:,:,2));
    %model(:,:,3) = 1 - sum(model,3);
end