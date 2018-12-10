function [] = fixTraits(traitMatrix,smoothValue)


    %% stack
    for e = 1:numel(traitMatrix)
        d = traitMatrix{e};
        for i = 1:size(d,1)
            for e2 = 1:(size(d,2))
                if strcmp(d{i,e2},'NA') | isempty(d{i,e2})
                    d{i,e2} = NaN;
                end
            end
        end
        d = cell2mat(d);
    end
    TM = d;
    %% remove NAN
    rm = any(isnan(TM),2);
    badTM = TM;
    TM(rm,:) = [];
    %% smooth
    TM = imfilter(TM,fspecial('average',[1 smoothValue]),'replicate');
    dTM = diff(TM,1,2);
    Z = zscore(dTM(:),1);
    Z = reshape(Z,size(dTM));
    rm = any(Z > 1,2)
    TM(any(rm,2),:) = [];
    %% perform PCA 
    [S C U E L ERR LAM] = PCA_FIT_FULL(TM,3);
    %% multple regressors
    
    
%%

    for e = 1:size(badTM,1)
        if sum(isnan(badTM(e,:))) < 5 && sum(isnan(badTM(e,:))) > 0
            rawTM = badTM(e,:);
          
            
            % PLS for NAN
            nidx = find(isnan(rawTM));
            nnidx = find(~isnan(rawTM));
            Y = TM(:,nidx);
            X = TM(:,nnidx);
            [XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(X,Y,4);
            fixedValues = [1 rawTM(nnidx)]*BETA;
            rawTMf = rawTM;
            rawTMf(nidx) = fixedValues;
            
            
            base = rawTMf;
            for sl = 1:1
                for l = 1:10
                    for r = 1:size(TM,2)
                        nr = setdiff(1:size(TM,2),r);
                        Y = TM(:,r);
                        X = TM(:,nr);
                        [XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(X,Y,3);
                        newCurve(r) = [1 rawTMf(nr)]*BETA;
                    end
                    dist(l) = norm(newCurve-base);
                    store(l,:) = newCurve;
                    rawTMf = newCurve;
                end
                [~,midx] = min(log(dist));
                newCurve = store(midx,:);


                %sim = newCurve;


                [xg(e,:),fval,exitflag,output] = fminunc(@(x)myCurveDistance(x,E,U,newCurve),mean(C,1));
                rawTMf = PCA_BKPROJ(xg(e,:),E,U);
            end
            sim = rawTMf;

            
            plot(rawTM,'k')
            hold on
            plot(sim,'r')
            %axis([0 13 0 1500])
            title(num2str(sum(isnan(badTM(i,:)))))
            waitforbuttonpress
            hold off
            drawnow
            
        end
    end
end