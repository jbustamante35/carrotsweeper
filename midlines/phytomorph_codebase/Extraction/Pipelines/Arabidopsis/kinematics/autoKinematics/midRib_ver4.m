function [midlineM] = midRib_ver4(FT,MSK,F1,F2)

    sFT = imfilter(FT,fspecial('average',[21 21]),'replicate');
    [S C U E L ERR LAM] = PCA_FIT_FULL(sFT,1);
    PV = 500;
    SCALE = 120;
    SIGMA = 90;
    STACK = zeros([size(sFT) numel(SCALE)]);
    ker = pdf('norm',linspace(-SCALE,SCALE,2*SCALE),0,SIGMA);
    ker = ones(1,300);
    for e = 1:size(sFT,1)
        sig = sFT(e,:);
        sig = bindVec(-sig);
        T = imfilter(sig,ker,'replicate');
        STACK(e,:,:) = permute(T,[2 1]);
    end
    
    f = mean(STACK,1);
    STACK = bsxfun(@times,STACK,f);
    
    
    %%% post process wavelet type 4
    midx = find(any(MSK,2));
    cnt = 1;
    mi = [];
    for e = 1:numel(midx)
        tmp = STACK(midx(e),:);
        [J,Midx] = max(tmp);
        mi(cnt,:) = [midx(e) Midx];
        cnt = cnt + 1;
    end
    midlineM = round(mi);
    midline = sub2ind([size(F1,1) size(F1,2)],midlineM(:,1),midlineM(:,2));
    midlineM = [];
    midlineM(:,1) = F1(midline);
    midlineM(:,2) = F2(midline);
    
    
end