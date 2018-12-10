function [R] = GO_SUM(S,holdOutGroups,type,display)
    
    R.compPreY = [];
    R.compAccY = [];
    
    R.compPreGrpMeanY = [];
    R.compAccGrpMeanY = [];
        
    R.dataPreY = [];
    R.dataAccY = [];
    
    R.dataPreGrpMeanY = [];
    R.dataAccGrpMeanY = [];

    R.classOut = {};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X = S.specData;
    Z = S.shapeData;
    Y = S.tipAngle;
    G = S.genoType;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform pca on X
    if type.PCAX.perform
        [xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,type.PCAX.dim);
    else
        xC = X;
        xU = mean(xC,1);
        xC = bsxfun(@minus,xC,xU);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform pca on Z
    if type.PCAZ.perform
        [zS zC zU zE zL zERR zLAM] = PCA_FIT_FULL(Z,type.PCAZ.dim);
    else
        zC = Z;
        zU = mean(zC,1);
        zC = bsxfun(@minus,zC,zU);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform pca on Y
    if type.PCAY.perform
        [yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,type.PCAY.dim);
    else
        yC = Y;,'LineWidth',5
        yU = mean(yC,1);
        yC = bsxfun(@minus,yC,yU);
    end
    
    % direct sum of vectors
    xC = [xC zC];
    
    
    UQ = unique(holdOutGroups);
    for u = 1:numel(UQ)        
        fidxTR = find(~strcmp(holdOutGroups,UQ{u}));
        fidxTE = find(strcmp(holdOutGroups,UQ{u}));
        
        TrainX = xC(fidxTR,:);
        TrainY = yC(fidxTR,:);
        
        
        if type.learnMethod.means
            TrainX = [];
            TrainY = [];
            genoGroups = S.genoType(fidxTR);
            UQg = unique(genoGroups);
            for ui = 1:numel(UQg)
                fidx = find(strcmp(UQg{ui},S.genoType));
                TrainX = [TrainX;mean(xC(fidx,:),1)];
                TrainY = [TrainY;mean(yC(fidx,:),1)];
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % learn
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch type.learnMethod.which
            case 'pls'            
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % simple pls regress
                [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(TrainX,TrainY,type.learnMethod.numComp);
                predictor = BETA(2:end,:);                
            case 'cca'
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % simple canoncorr
                [A,B,r,U,V,stats] = canoncorr(TrainX,TrainY);
                [mA mB] = myCCA(xC(fidxTR,:),yC(fidxTR,:),3);
                predictor = A*inv(B);
            case 'net'
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % neural network
                type.learnMethod.net = feedforwardnet([30]);
                [type.learnMethod.net tr] = train(type.learnMethod.net,TrainX',TrainY');
            case 'kMani'
                lM = lManifold();
                lM.setmodelCompX(type.PCAX.dim);
                lM.setmodelCompY(type.PCAY.dim);
                lM.addXY(TrainX,TrainY);
                lM.setGroupN(3);
                lM.learn();
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store which group
        for e = 1:numel(fidxTE)
            R.classOut{end+1} = UQ{u};
        end
        R.classOutOrder = UQ;
        
        switch type.learnMethod.which            
            case 'pls'
                % perform prediction
                prediction = xC(fidxTE,:)*predictor;                
            case 'cca'
                % perform prediction
                prediction = xC(fidxTE,:)*predictor;
            case 'net'
                prediction = type.learnMethod.net(xC(fidxTE,:)')';
            case 'kMani'
                prediction = lM.predict(xC(fidxTE,:));
        end
        
        % if pca the store comps and backproject into larger space
        if type.PCAY.perform
            % store prediction on comps
            R.compPreY = [R.compPreY;prediction];
            R.compAccY = [R.compAccY;yC(fidxTE,:)];
            % store mean comps for holdout group
            R.compPreGrpMeanY = [R.compPreGrpMeanY;mean(prediction,1)];
            R.compAccGrpMeanY = [R.compAccGrpMeanY;mean(yC(fidxTE,:),1)];
            % back project
            prediction = PCA_BKPROJ(prediction,yE,yU);
        end
        
        % store data predictions
        R.dataPreY = [R.dataPreY;prediction];
        R.dataAccY = [R.dataAccY;Y(fidxTE,:)];
        
        % store data predictions
        R.dataPreGrpMeanY = [R.dataPreGrpMeanY;mean(prediction,1)];
        R.dataAccGrpMeanY = [R.dataAccGrpMeanY;mean(Y(fidxTE,:),1)];
        
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % display
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % look at pca predictions 
    if display.pca.perform
        for e = 1:size(R.compPreY,2)
            subX = R.compPreY(:,e);
            subY = R.compAccY(:,e);
            [RHO PVAL] = corr(subX,subY);
            MX = max(([subX(:);subY(:)]));
            MN = min(([subX(:);subY(:)]));
            figure;
            plot(subX,subY,'.');            
            hold on
            plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
            ylabel('Predicted');
            xlabel('Actual');
            axis([MN MX MN MX]);
            title([display.pca.domain '-->' display.pca.codomain '- PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
            %saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % look at pca predictions 
    
    if display.pca.perform
        for e = 1:size(R.compPreY,2)
            subX = R.compPreGrpMeanY(:,e);
            subY = R.compAccGrpMeanY(:,e);
            [RHO PVAL] = corr(subX,subY);
            MX = max(([subX(:);subY(:)]));
            MN = min(([subX(:);subY(:)]));
            figure;
            for i = 1:size(subX,1)
                plot(subX(i),subY(i),'.');
                hold all
            end
            axis([MN MX MN MX]);            
            plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
            ylabel('Predicted');
            xlabel('Actual');
            title([display.pca.domain '-->' display.pca.codomain '- PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);            
            %saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
            legend(R.classOutOrder);
        end
    end
    
    UQg = unique(R.classOut);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % look at data predictions 
    if display.indivdual.perform
        scale = 180/pi;
        figure;
        for e = 1:numel(UQg)
            fidx = find(strcmp(UQg{e},R.classOut));
            if fidx > 1
                preY = R.dataPreY(fidx,:);
                accY = R.dataAccY(fidx,:);
                uP = mean(preY,1);
                sP = std(preY,1,1);
                uA = mean(accY,1);
                sA = std(accY,1,1);
                plot(scale*preY','m--','LineWidth',1);
                hold on
                plot(scale*accY','k--','LineWidth',1);
                errorbar(scale*uP,scale*sP,'r','LineWidth',3);
                hold on
                errorbar(scale*uA,scale*sA,'k','LineWidth',3);                
                hold off
                axis([0 61 -10 100]);
                title(UQg{e});
                drawnow
                pause(1);
                %title([display.pca.domain '-->' display.pca.codomain '- PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
                %saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
            end
        end
    end
    
end