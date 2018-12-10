function [f,phase,fIDX,sigRAW] = findMajorFrequency(subI,cutoff,fCut,fIDX,disp,F,FF)


    toPROC = bsxfun(@minus,subI,mean(subI,2));
    BK = imfilter(mean(toPROC,1),fspecial('average',[1 400]),'replicate');
    toPROC = bsxfun(@minus,subI,BK);
    toPROC = bsxfun(@minus,subI,mean(subI,2));
    sig = fft(toPROC,[],2);
    sig2 = mean(abs(sig),1);
    
    
    %{
    if ~isempty(F)
        sig2 = -log(normpdf(sig2,F,FF));
    end
    %}
    
    sigRAW = sig2;
    ssig2 = imfilter(sig2,fspecial('average',[1 41]),'replicate');
    if disp
        figure;
        plot(ssig2,'k')
        hold on
        plot(sig2);
        figure;
    end
    
    sig3 = sig2 - ssig2;
    if disp
        plot(sig3);
    end
    % change dilate amount from 51 to 31
    cutoff = round(size(subI,2)*fCut/2);
    sig3(1:cutoff) = 0;
    sig4 = imdilate(sig3,strel('disk',31,0)) == sig3;
    
    
    
    if isempty(fIDX)
        fidx = find(sig4);
    else
        fidx(1) = fIDX;
    end
    if disp
        hold on
        plot(sig4*30);
    end
    %fidx = fidx - 1;
    if ~isempty(fidx)
        f = (size(subI,2)/(fidx(1)-1))^-1;
        phase = angle(sig(:,fidx(1)));
    else
        f = 0;
        phase = zeros(size(subI,1),1);
    end
    
    
    
    
end

%{
    f = findMajorFrequency(subI,50);

%}