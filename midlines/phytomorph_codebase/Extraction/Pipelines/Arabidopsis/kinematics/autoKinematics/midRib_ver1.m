function [midlineM] = midRib_ver1(FT,MSK,F1,F2)

    sFT = imfilter(FT,fspecial('average',[21 21]),'replicate');
    %[S C U E L ERR LAM] = PCA_FIT_FULL(sFT,1);
    PV = 500;
    SCALE = 120;
    STACK = zeros([size(sFT) numel(SCALE)]);
    for e = 1:size(sFT,1)
        sig = sFT(e,:);
        sig = padarray(sig,[0 PV],'replicate');
        T = cwt(sig,SCALE,'gaus1');
        T(:,1:PV) = [];
        T(:,(end-(PV-1)):end) = [];
        STACK(e,:,:) = permute(T,[2 1]);
    end
    
    %%% post process wavelet type 2
    midx = find(any(MSK,2));
    cnt = 1;
    mi = [];
    for e = 1:numel(midx)
        tmp = STACK(midx(e),:);
        fidx_max = find(imdilate(tmp,ones([1 20])) == tmp);
        fidx_min = find(imdilate(-tmp,ones([1 20])) == -tmp);
        [J,m1] = max(tmp(fidx_max));
        [J,m2] = max(-tmp(fidx_min));
        fidx_max = fidx_max(m1);
        fidx_min = fidx_min(m2);
        mi(cnt,:) = [midx(e) mean([fidx_max fidx_min])];
        cnt = cnt + 1;
    end
    midlineM = round(mi);
    midline = sub2ind([size(F1,1) size(F1,2)],midlineM(:,1),midlineM(:,2));
    midlineM = [];
    midlineM(:,1) = F1(midline);
    midlineM(:,2) = F2(midline);
    
    
end