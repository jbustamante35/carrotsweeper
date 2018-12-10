function [beta,C] = findRegressor(rX,rY,L)
    disp = 0;
    for l = 1:L
        parfor e = 1:size(rX,1)
            idx = setdiff(1:size(rX,1),e);
            % subset the data
            subX = rX(idx,:);
            subY = rY(idx);
            % fit function
            [Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse,stats,Weights] = plsregress(subX,subY,l);
            % predict
            Yp(l,e) = ([1 rX(e,:)]*beta);
        end
        % correlate over prediction
        COR(l) = corr(rY(:),Yp(l,:)');
        if disp
            figure(h1)
            plot(COR)
            drawnow
        end
    end
    
    [C,loc] = max(real(COR));
    [Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse,stats,Weights] = plsregress(rX,rY,loc);
    
end