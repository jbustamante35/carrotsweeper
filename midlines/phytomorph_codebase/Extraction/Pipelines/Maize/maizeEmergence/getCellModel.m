function [idx,model,out] = getCellModel(I,obj,disp)
    sz = size(I);
    zMASK = I(:,:,1) ~= 0;
    out = I;
    I = reshape(I,[prod(sz(1:2)) sz(3)]);
    [idx,~,prob] = cluster(obj,I);
    idx = reshape(idx,sz(1:2));
    idx = idx.*double(zMASK);
    CL = {'b' 'r' 'b'};
    
    for k = 1:3
        out(:,:,k) = bindVec(out(:,:,k));
    end
    for k = 1:size(prob,2)
        model(:,:,k) = reshape(prob(:,k),sz(1:2));
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