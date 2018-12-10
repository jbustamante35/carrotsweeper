function [IDX M] = matchScale(I,vec,patchSize,reSize,E,U)
    for r = 1:numel(reSize)
        oI = imresize(I,reSize(r));
        rI = oI - imfilter(oI,fspecial('average',patchSize));
        B = [];
        for k = 1:size(rI,3)
            B = [B ; im2colF(rI(:,:,k),patchSize,[1 1])];
        end
        if nargin > 4
            B = PCA_REPROJ(B',E,U)';
        end
        
        
        B = [ones(1,size(B,2)) ;B];
        img = vec'*B;
        H = col2im(img,patchSize,size(rI(:,:,1)));
        
        [dx dy] = find(H == imdilate(H,strel('disk',31,0)));
        [idx] = find(H == imdilate(H,strel('disk',31,0)));
        [~,sidx] = sort(H(idx),'descend');
        dx = dx(sidx);
        dy = dy(sidx);
        [M(r),midx] = max(H(idx));
        IDX(r,:) = [dy(midx) dx(midx)];
        imshow([cat(3,bindVec(H),bindVec(H),bindVec(H)),oI((patchSize(1)-1)/2:(end-(patchSize(1)-1)/2-1),(patchSize(2)-1)/2:(end-(patchSize(2)-1)/2-1),:)],[]);
        hold on
        plot(dy(1:min(5,numel(sidx)))+size(H,2),dx(1:min(5,numel(sidx))),'mo');
        drawnow
        hold off
    end
end