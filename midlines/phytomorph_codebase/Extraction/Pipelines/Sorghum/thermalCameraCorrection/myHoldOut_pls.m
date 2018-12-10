function [CORR] = myHoldOut_pls(rX,rY,L,R,holdOutMethod,M)
    h1 = figure;
    h2 = figure;
    
    COR = [];
    Yp = [];
    
    % sweep over the L models
    for l = 1:L
        
        tmpStore = zeros(1,size(rY,2));
        
        parfor e = 1:R
            
            [Train, Test] = crossvalind(holdOutMethod,size(rX,1),M);
            
            

            subX = rX(Train,:);
            subY = rY(Train,:);


            
            [Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse,stats,Weights] = plsregress(subX,subY,l);
            Yp = [ones(sum(Test),1) rX(Test,:)]*beta;
            tmpStore = tmpStore + diag(corr(Yp,rY(Test,:)))';
        end
        
        tmpStore = tmpStore/R;
        CORR(l,:) = tmpStore;
    end
end