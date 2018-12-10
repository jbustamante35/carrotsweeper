function [SAM] = gridSample(I,curveBundle,frameBundle,distance,nP)    
    SAM = zeros([nP 3]);
    for m = 1:nP(2)
        SP = curveBundle(m,:) - distance*frameBundle(m,:,2);
        EP = curveBundle(m,:) + distance*frameBundle(m,:,2);
        lSEG = [linspace(SP(1),EP(1),nP(1));linspace(SP(2),EP(2),nP(1))]';
        SAM(:,m,1) = ba_interp2(I,lSEG(:,2),lSEG(:,1));
        SAM(:,m,2:3) = lSEG;
    end            
end