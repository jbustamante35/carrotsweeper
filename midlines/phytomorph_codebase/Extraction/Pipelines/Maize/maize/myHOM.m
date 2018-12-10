function [CORR masterX masterY mU mV RMSEP RMSEC] = myHOM(X,Y,labels,nx,ny)
    UQ = unique(labels);
    [xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,nx);
    [yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,ny);
    for u = 1:numel(UQ)        
        subX = xC;
        subY = yC;
        trxC = subX(strcmp(UQ,UQ{u}),:);
        tryC = subY(strcmp(UQ,UQ{u}),:);
        subX(strcmp(UQ,UQ{u}),:) = [];
        subY(strcmp(UQ,UQ{u}),:) = [];
        xU = mean(subX);
        yU = mean(subY);
        [mA,mB,mr,mU,mV,mstats] = canoncorr(subX,subY);
        masterX(u,:) = (trxC-xU)*mA;
        masterY(u,:) = (tryC-yU)*mB;
        
        
        
        predictY = (inv(mB')*mU')';
        delta = (predictY - subY);
        RMSEC(u) = mean(sum(delta.*delta,2).^.5);
        
        predictY = (inv(mB')*masterX(u,:)')';
        delta = (predictY - tryC);
        RMSEP(u) = mean(sum(delta.*delta,2).^.5);
        
    end
    CORR = corr(masterX,masterY);
    
    
    
    
    
    
    [mA,mB,mr,mU,mV,mstats] = canoncorr(xC,yC);
    
    
    
    
    
    
    %{
    for e = 1:size(mU,2)
        figure;
        plot(mU(:,e),mV(:,e),'.');
        title(num2str(corr(mU(:,e),mV(:,e))));            
        hold on
        plot(linspace(-4,4,2),linspace(-4,4,2))
    end
    %}
end