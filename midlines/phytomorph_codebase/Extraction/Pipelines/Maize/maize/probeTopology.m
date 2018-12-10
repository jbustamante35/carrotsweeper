function [EXP PVAL tVecMax] = probeTopology(networkFunction,maxNeurons,maxLayers,tr,X,Y,Labels,E,U)
    UQ = unique(Labels);
    for layerNumber = 1:maxLayers
        tVecMax = maxNeurons*ones(1,layerNumber);        
        [net TOPO] = generateNetwork(tVecMax,networkFunction);
        for t = 1:tr
            r = [];p = [];
            for u = 1:numel(UQ)
                net = networkFunction(TOPO);
                %net.trainParam.showWindow = false;
                fidx = ~strcmp(Labels,UQ{u});
                tidx = strcmp(Labels,UQ{u});
                subX = X(fidx,:);
                subY = Y(fidx,:);
                testX = X(tidx,:);
                testY = Y(tidx,:);
                [RAW PRE] = testNetwork(subX,subY,testX,testY,E,U,net);
                r = [r;RAW];
                p = [p;PRE];
            end            
            [COR PAL] = corr(p,r);
            EXP(t,:) = diag(COR);
            PVAL(t,:) = diag(PAL);
        end
    end
end

function [prediction netPrediction] = testNetwork(XTrain,YTrain,XTest,YTest,E,U,net)
    net = train(net,{XTrain'},{YTrain'});
    Y = net(XTest');
    prediction = PCA_BKPROJ(YTest,E,U);
    netPrediction = PCA_BKPROJ(Y',E,U);
    %{
    errorbar(mean(prediction,1),std(prediction,1,1),'k');
    hold on
    errorbar(mean(netPrediction,1),std(netPrediction,1,1),'r');
    hold off
    drawnow
    %}
end


function [net TOPO] = generateNetwork(tVecMax,networkFunction)
    TOPO = [];
    for e = 1:numel(tVecMax)
        TOPO = [TOPO randi(tVecMax(e),1)];
    end
    net = networkFunction(TOPO);
end