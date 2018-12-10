function [midlineM] = midRib_ver0(FT,MSK,F1,F2)

    sFT = imfilter(FT,fspecial('average',[21 21]),'replicate');
    %[S C U E L ERR LAM] = PCA_FIT_FULL(sFT,1);
    PV = 500;
    SCALE = 60;
    STACK = zeros([size(sFT) numel(SCALE)]);
    for e = 1:size(sFT,1)
        sig = sFT(e,:);
        sig = padarray(sig,[0 PV],'replicate');
        T = cwt(sig,SCALE,'gaus1');
        T(:,1:PV) = [];
        T(:,(end-(PV-1)):end) = [];
        STACK(e,:,:) = permute(T,[2 1]);
    end
    
    %%% post process wavelet type 1
    PDF = zeros(size(FT));
    for e = 1:size(STACK,3)
        sig = squeeze(STACK(:,:,e));
        sig = abs(sig);
        thresh = graythresh(sig);
        tmpM = sig > thresh;
        PDF = PDF + tmpM;
    end
    PDF = imfilter(PDF,fspecial('average',[21 21]),'replicate');
    tmpM = PDF > graythresh(PDF);
    tmpM = tmpM.*MSK;
    tmpV = tmpM(end,:);
    tmpM(end,:) = 1;
    %tmpM = bwlarge(tmpM);
    tmpM(end,:) = tmpV;
    midlineF = [];
    for e = 1:size(tmpM,1)
        sig = tmpM(e,:);
        fidx = find(sig);
        if ~isempty(fidx)
            midlineF = [midlineF;[e mean([fidx(1) fidx(end)])]];
        end
    end
    SNIP = 10;
    midlineF(1:SNIP,:) = [];
    midlineF = imfilter(midlineF',fspecial('average',[1 100]),'replicate');
    midlineF = arcLength(midlineF','arcLen')';
    midlineF = midlineF';
    midlineM = midlineF;
    midlineM = round(midlineM);
    midline = sub2ind([size(F1,1) size(F1,2)],midlineM(:,1),midlineM(:,2));
    midlineM = [];
    midlineM(:,1) = F1(midline);
    midlineM(:,2) = F2(midline);
    
end