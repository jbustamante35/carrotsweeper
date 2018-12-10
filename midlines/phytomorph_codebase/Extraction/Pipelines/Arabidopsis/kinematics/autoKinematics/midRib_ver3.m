function [midlineM] = midRib_ver3(FT,MSK,F1,F2)
    sFT = imfilter(FT,fspecial('disk',[71]),'replicate');
    REV = 300-1;
    dist = [];
    sFT = padarray(sFT,[0 REV],'replicate','post');
    dist = zeros(size(FT));
    for e = 1:size(sFT,1)
        sig = sFT(e,:);
        sig = im2col(sig,[1 REV+1]);
        UPPER = sig(1:(end/2),:);
        LOWER = flipud(sig((end/2+1):end,:));
        d = sum((UPPER - LOWER).^2,1).^.5;
        dist(e,:) = d';
        fprintf(['Done with row:' num2str(e) ':' num2str(size(sFT,1)) '\n']);
    end
    
    dist = bsxfun(@times,dist,mean(dist,1));
    %dist = dist.*MSK;
    %%% post process wavelet type 4
    widx = find(any(MSK,1))-(REV+1)/2;
    SNIP = round(numel(widx)/4);
    widx = widx(SNIP:(end-SNIP));
    cnt = 1;
    mi = [];
    midx = find(any(MSK,2));
    for e = 1:numel(midx)
        tmp = dist(midx(e),widx);
        %plot(tmp)
        %drawnow
        [J,Midx] = min(tmp);
        mi(cnt,:) = [midx(e) Midx+widx(1)+(REV+1)/2];
        cnt = cnt + 1;
       % waitforbuttonpress
    end
    midlineM = round(mi);
    midline = sub2ind([size(F1,1) size(F1,2)],midlineM(:,1),midlineM(:,2));
    midlineM = [];
    midlineM(:,1) = F1(midline);
    midlineM(:,2) = F2(midline);
    
    
end