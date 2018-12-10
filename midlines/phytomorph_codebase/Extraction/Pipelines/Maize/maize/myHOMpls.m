function [PRESS] = myHOMpls(X,Y,labels,factors)
    

    
    [X,muX,sigmaX] = zscore(X);
    [Y,muY,sigmaY] = zscore(Y);
    %muY = mean(Y,1);
    %X = bsxfun(@plus,X,-mean(X,1));
    %Y = bsxfun(@plus,Y,-mean(Y,1));
    
    UQ = unique(labels);
    %[xS X xU xE xL xERR xLAM] = PCA_FIT_FULL(X,61);
    %[yS Y yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,61);
    
    
    
    UQ = unique(labels);
    for u = 1:numel(UQ)
        tic
        subX = X;
        subY = Y;
        trxC = subX(strcmp(UQ,UQ{u}),:);
        tryC = subY(strcmp(UQ,UQ{u}),:);
        subX(strcmp(UQ,UQ{u}),:) = [];
        subY(strcmp(UQ,UQ{u}),:) = [];
        
        options = statset('UseParallel','always');
        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(subX,subY,factors,'Options',options,'cv',10);
        %[wx,wy,sx,sy,PRE] = myPLS_EZ(subX,subY,factors);
        %Ypredict = (wy*diag(PRE)*(trxC*wx)')';
        Ypredict = [ones(1,1),trxC]*BETA;
        %Ypredict = bsxfun(@times,Ypredict,sigmaY.^-1);
        %Ypredict = bsxfun(@plus,Ypredict,muY);
        %Yacc = bsxfun(@plus,tryC,muY);
        PRESS(u) = norm(Ypredict-tryC)^2;
        
        %Ypredict = PCA_BKPROJ(Ypredict,yE,yU);
        %Yacc = PCA_BKPROJ(tryC,yE,yU);
        %{
        plot(Ypredict,'r')
        hold on
        plot(Yacc,'b')
        drawnow
        pause(.2);
        hold off
        %}
        
        %{
        masterX(u,:) = (trxC-xU)*mA;
        masterY(u,:) = (tryC-yU)*mB;
        
        
        
        predictY = (inv(mB')*mU')';
        delta = (predictY - subY);
        RMSEC(u) = mean(sum(delta.*delta,2).^.5);
        
        predictY = (inv(mB')*masterX(u,:)')';
        delta = (predictY - tryC);
        RMSEP(u) = mean(sum(delta.*delta,2).^.5);
        %}
        u
        numel(UQ)
        toc
    end
    
    
    
    
    
    
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