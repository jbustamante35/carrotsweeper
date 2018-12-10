function [vB1,hB1,vB2,hB2,simBK] = generateBKsignals(I,BK)
    [vB1,hB1] = h0(I,BK);
    [hB2,vB2] = h0(permute(I,[2 1 3]),permute(BK,[2 1 3]));
    vB1 = shiftdim(vB1',-1);
    hB1 = reshape(hB1,[size(hB1,1) 1 3]);
    
    vB2 = vB2';
    vB2 = shiftdim(vB2',-1);
    hB2 = reshape(hB2',[size(hB2,2) 1 3]);
    
    [simBK] = generateBackground(vB1,hB1,vB2,hB2,size(I));
    
end

function [hB,vB] = h0(I,BK)
        mI = bsxfun(@times,I,BK);
        hB = [];
        hSZ = 401;
        for k = 1:size(mI,3)
            tmp = mI(:,:,k);
            tmp(tmp==0) = NaN;
            tmp = nanmean(tmp,1);
            tmp(isnan(tmp)) = nanmean(tmp);
            hB(k,:) = tmp;
            hB(k,:) = imfilter(hB(k,:),fspecial('average',[1 hSZ]),'replicate');
        end

        tmpBK = shiftdim(hB',-1);
        tmpBK = repmat(tmpBK,[size(I,1) 1 1]);

        deltaBK = I - tmpBK;

        deltaBK = bsxfun(@times,deltaBK,BK);
        vB = [];
        vSZ = 401;
        for k = 1:size(deltaBK,3)
            tmp = deltaBK(:,:,k);
            tmp(tmp==0) = NaN;
            tmp = nanmean(tmp,2);
            tmp(isnan(tmp)) = nanmean(tmp);
            vB(:,k) = tmp;
            vB(:,k) = imfilter(vB(:,k) ,fspecial('average',[vSZ 1]),'replicate');
        end
        DtmpBK = shiftdim(vB,-1);
        DtmpBK = repmat(DtmpBK,[size(I,2) 1 1]);
        DtmpBK = permute(DtmpBK,[2 1 3]);
end