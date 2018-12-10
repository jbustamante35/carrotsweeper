function [] = findExpRegressorFunction(rX,rY,n,disp)

    rX = log(rX);
    for l = 1:10
        parfor e = 1:size(rX,1)
            idx = setdiff(1:size(rX,1),e);
            % subset the data
            subX = rX(idx,:);
            subY = rY(idx);
            % fit function
            [Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse,stats,Weights] = plsregress(subX,subY,l);
            % predict
            Yp(l,e) = [1 rX(e,:)]*beta;
        end
        % correlate over prediction
        COR(l) = corr(rY(:),Yp(l,:)','type','Pearson');
        if disp
            figure(h1)
            plot(COR)
            drawnow
        end
    end
end


