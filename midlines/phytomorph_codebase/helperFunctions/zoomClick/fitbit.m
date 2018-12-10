function [map,beta,WINDOW] = fitbit(X,Y,bx,RNG,map)

    if nargin == 4
        map = [];
    end
    
    
    [dX,map] = reduce(X,bx,map);
    %RGB = ind2rgb(dX(:,:,2^bx),map);
    [bX] = bits(dX,bx);
    
    
    
    bX = double(bX);
    B = im2colF((bX(:,:,1)),[15 15],[1 1]);
    B = zeros(size(B,1),size(B,2),size(dX,3));
    for tr = 1:size(bX,3)
        %tic;
        B(:,:,tr) = im2colF((bX(:,:,tr)),[15 15],[1 1]);
        %tm(tr) = toc;
        %mean(tr)*(size(dX,3)-tr)/60
        tr
    end
    B = squeeze(sum(B,1));
    bX = B';
    
    for r = 1:size(Y,2)
        WINDOW(:,r) = linspace(min(Y(:,r),[],1),max(Y(:,r),[],1),RNG)';
        for level = 1:size(WINDOW,1)
            pY(:,level,r) = Y(:,r) > WINDOW(level,r);
        end
    end
    
    
    for r = 1:size(Y,2)
        [Xloadings,Yloadings,Xscores,Yscores,beta(:,:,r),pctVar,mse,stats,Weights] = plsregress(bX,pY(:,:,r),11);
        r
    end

end