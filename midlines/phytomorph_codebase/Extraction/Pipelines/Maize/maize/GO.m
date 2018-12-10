function [R] = GO(S,holdOutGroups,type,display)
    

    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % assign X and Y
    X = S.specData;    
    %p = randperm(size(X,1));
    %X = X(p,:);
    Y = S.tipAngle;
    G = S.genoType;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pre-process
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    % perform pca on Y
    if type.PCAY.perform
        [yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,type.PCAY.dim);
    else
        yC = Y;
        yU = mean(yC,1);
        yC = bsxfun(@minus,yC,yU);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % learn and predict for each level of each grouping factor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for f = 1:numel(S.groupingFactors)
        if S.groupingFactors(f).toGroup
            % get each level for the fth grouping factor
            UQ_levels = unique(S.groupingFactors(f).Groups);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % learn and predict -START
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
            for u = 1:numel(UQ_levels)
                % get the index for the training set
                fidxTR = find(~strcmp(S.groupingFactors(f).Groups,UQ_levels{u}));
                % get the index for the testing set
                fidxTE = find(strcmp(S.groupingFactors(f).Groups,UQ_levels{u}));
                fidxTE = 1:numel(S.groupingFactors(f).Groups);
                % get the train setX
                TrainX = xC(fidxTR,:);
                % get the train setY
                TrainY = yC(fidxTR,:);

                %{
                % to to learn on factor means
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
                %}


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % learn - START
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
                        type.learnMethod.net = feedforwardnet([15 5 15]);
                        %type.learnMethod.net = newrb(TrainX',TrainY',0,10,round(size(TrainX,1)),50);
                        [type.learnMethod.net tr] = train(type.learnMethod.net,TrainX',TrainY');
                    case 'kMani'
                        lM = lManifold();
                        lM.setmodelCompX(type.PCAX.dim);
                        lM.setmodelCompY(type.PCAY.dim);
                        lM.addXY(TrainX,TrainY);
                        lM.setGroupN(3);
                        lM.learn();
                    case 'lookup'
                        learnSETX = TrainX;
                        learnSETY = TrainY;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % learn STOP
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


                
                %{
                % store which group
                for e = 1:numel(fidxTE)
                    for sf = 1:numel(S.groupingFactors)
                        R{f}{u}.groupingFactors(sf).Groups = S.groupingFactors(sf).Groups(fidxTE);
                    end
                end
                %}
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % predict START
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
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
                    case 'lookup'
                        prediction = [];
                        for e = 1:numel(fidxTE)
                            tmp = xC(fidxTE(e),:);
                            %tmp = tmp/norm(tmp);
                            %tmp = learnSETX*tmp';
                            tmp = bsxfun(@minus,learnSETX,tmp);
                            %tmp = sum(tmp.*tmp,2);
                            %tmp = sum(abs(tmp)*diag(xLAM)',2);
                            tmp = abs(tmp)*diag(xLAM).^-.5';
                            %tmp = sum(tmp.*tmp,2).^.5;
                            [tmp idx] = sort(tmp);
                            %idx = randi(size(learnSETY,1),1);
                            %W = tmp(1:3)/sum(tmp(1:3));
                            %prediction(e,:) = W'*learnSETY(idx(1:3),:);
                            %prediction(e,:) = mean(learnSETY(idx(1),:),1);
                            prediction(e,:) = mean(learnSETY(idx(1),:),1);
                        end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % predict
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % store results - START
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                R{f}{u}.compPreY = [];
                R{f}{u}.compAccY = [];

                R{f}{u}.compPreGrpMeanY = [];
                R{f}{u}.compAccGrpMeanY = [];

                R{f}{u}.dataPreY = [];
                R{f}{u}.dataAccY = [];

                R{f}{u}.dataPreGrpMeanY = [];
                R{f}{u}.dataAccGrpMeanY = [];

                R{f}{u}.classOut = {};
                
                % if pca the store comps and backproject into larger space
                if type.PCAY.perform
                    % store prediction on comps
                    R{f}{u}.compPreY = [R{f}{u}.compPreY;prediction];
                    R{f}{u}.compAccY = [R{f}{u}.compAccY;yC(fidxTE,:)];
                    % store mean comps for holdout group
                    R{f}{u}.compPreGrpMeanY = [R{f}{u}.compPreGrpMeanY;mean(prediction,1)];
                    R{f}{u}.compAccGrpMeanY = [R{f}{u}.compAccGrpMeanY;mean(yC(fidxTE,:),1)];
                    % back project
                    prediction = PCA_BKPROJ(prediction,yE,yU);
                end

                % store data predictions
                R{f}{u}.dataPreY = [R{f}{u}.dataPreY;prediction];
                R{f}{u}.dataAccY = [R{f}{u}.dataAccY;Y(fidxTE,:)];

                % store data predictions
                R{f}{u}.dataPreGrpMeanY = [R{f}{u}.dataPreGrpMeanY;mean(prediction,1)];
                R{f}{u}.dataAccGrpMeanY = [R{f}{u}.dataAccGrpMeanY;mean(Y(fidxTE,:),1)];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % store results - STOP
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % display
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for f = 1:numel(S.groupingFactors)
        if S.groupingFactors(f).toGroup
            for df = 1:numel(S.groupingFactors)
                if S.groupingFactors(df).toDisplay
                    % get each level for the fth grouping factor
                    UQ_levels = unique(S.groupingFactors(f).Groups);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % loop over each predicting factor
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    AccrosslevelName = {};
                    dataAccrosLevelsAcc = [];
                    dataAccrosLevelsPre = [];
                    for u = 1:numel(UQ_levels)
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % loop over each display factor
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % get each level for the fth grouping factor
                        dUQ_levels = unique(S.groupingFactors(df).Groups);
                        for du = 1:numel(dUQ_levels)
                            % if the prediction level == the display level
                            if strcmp(UQ_levels{u},dUQ_levels{du})
                                idx = find(strcmp(dUQ_levels{du},S.groupingFactors(df).Groups));
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % look at pca predictions - indivduals
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                if display.pca.perform.level | display.pca.perform.across
                                    pY = R{f}{u}.compPreY;
                                    pX = R{f}{u}.compAccY;
                                    % loop over princ comps
                                    for e = 1:size(R{f}{u}.compPreY,2)
                                        subX = pX(idx,e);
                                        subY = pY(idx,e);
                                        [RHO PVAL] = corr(subX,subY);
                                        MX = max(([subX(:);subY(:)]));
                                        MN = min(([subX(:);subY(:)]));
                                        if display.pca.perform.level
                                            figure;
                                            plot(subX,subY,'.');            
                                            hold on
                                            plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
                                            ylabel('Predicted');
                                            xlabel('Actual');
                                            axis([MN MX MN MX]);
                                            title([{['Factor Predict->' S.groupingFactors(f).Name ...
                                                   ' Level Predict->' UQ_levels{u}]}, ...
                                                   {['Factor Display->' S.groupingFactors(df).Name ...
                                                   ' Level Display->' dUQ_levels{du}]}, ...
                                                   {[display.pca.domain '-->' display.pca.codomain '- PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]}]);
                                            %saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
                                        end
                                    end
                                    % store for display accross factors
                                    for i = 1:numel(idx)
                                        AccrosslevelName{end+1} = dUQ_levels{du};
                                    end
                                    dataAccrosLevelsAcc = [dataAccrosLevelsAcc;pX];
                                    dataAccrosLevelsPre = [dataAccrosLevelsPre;pY];
                                end
                            end
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % display across level for factor
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % look at predictions via the groupings
                    if display.pca.perform.across
                        CL = {'c*','r*','k*','b*','y*','m*','c^','r^','k^','b^','y^','m^','c.','r.','k.','b.','y.','m.','co','ro','ko','bo','yo','mo','cs','rs','ks','bs','ys','ms'};
                        UQL = unique(AccrosslevelName);               
                        for ee = 1:size(dataAccrosLevelsAcc,2)
                            figure;
                            hold all;
                            % for each group in the display grouping
                            for du = 1:numel(UQL)
                                idx = strcmp(AccrosslevelName,UQL{du});
                                subX = dataAccrosLevelsAcc(idx,ee);
                                subY = dataAccrosLevelsPre(idx,ee);
                                plot(mean(subY,1),mean(subX,1),CL{du},'MarkerSize',10,'MarkerFaceColor',CL{du}(1));
                            end
                            legend(UQL);
                            % for each group in the display grouping
                            for du = 1:numel(UQL)
                                idx = strcmp(AccrosslevelName,UQL{du});
                                subX = dataAccrosLevelsAcc(idx,ee);
                                subY = dataAccrosLevelsPre(idx,ee);
                                plot(subY,subX,CL{du});
                            end
                            subX = dataAccrosLevelsAcc(:,ee);
                            subY = dataAccrosLevelsPre(:,ee);
                            [RHO PVAL] = corr(subX,subY);
                            MX = max(([subX(:);subY(:)]));
                            MN = min(([subX(:);subY(:)]));
                            axis([MN MX MN MX]);            
                            plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
                            ylabel('Predicted');
                            xlabel('Actual');
                            title([{['Factor Predict->' S.groupingFactors(f).Name ...
                                   ' Level Predict->' UQ_levels{u}]}, ...
                                   {['Accorss Factor Display->' S.groupingFactors(df).Name]}, ...
                                    {[ display.pca.domain '-->' display.pca.codomain '- PC' num2str(ee) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]}]);       
                        end
                    end
                    %{
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % look at data predictions 
                    if display.indivdual.perform
                        figure;
                        scale = 180/pi;
                        preY = R{f}{u}.dataPreY;
                        accY = R{f}{u}.dataAccY;
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
                        drawnow
                        pause(.3);
                        [RHO PVAL] = corr(preY,accY);
                        title([{S.groupingFactors(f).Name}, ...
                                   {['Level Predict->' UQ_levels{u}]}, ...
                                   {[display.pca.domain '-->' display.pca.codomain]}]);
                         %saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
                    end
                    %}
                end
            end
        end
    end
end