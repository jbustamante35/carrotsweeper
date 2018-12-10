function [tester] = shapeVerticalStripNozzle(imageNozzle,maskNozzle,nCOMP,transPROB,repeatingUnits,nGRPs,LOOP)
    try
        
        
        
        % collect full Y spray for the init Y
        Y = maskNozzle.collectFullSpray();
        nCOMP_in = nCOMP;
        nCOMP = 1;
        for loop = 1:LOOP
            
        


            % other X reduction techniques
            %reducedImageNozzle1 = imageNozzle.getTransductionNozzle(@(X)getPLSDR_func(X,Y,nCOMP));
            %reducedImageNozzle2 = imageNozzle.getTransductionNozzle(@(X)getLNAn_func(X,Y,nCOMP));
            %reducedImageNozzle3 = imageNozzle.getReductionNozzle(nCOMP);
            %reducedImageNozzle = imageNozzle.getTransductionNozzle(@(X)getLNAn_func(X,Y,nCOMP));
            %reducedImageNozzle = imageNozzle.getTransductionNozzle(@(X)getPLSDR_func(X,Y,nCOMP));


            reducedImageNozzle = imageNozzle.getTransductionNozzle(@(X)getPNN_func(X,Y,nCOMP_in));


            % reduce the dim of the spray
            %reducedImageNozzle = imageNozzle.getReductionNozzle(nCOMP);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % hard code to two classes
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for e = 1:2
                % plug nozzles into nozzleManifold to create filter-nozzle
                filteredClassNozzle{e} = nozzleManifold({reducedImageNozzle,maskNozzle},@(X,e0,e1)filterSpray(X,e-1),1);
                % hook up mean gauge to filter nozzle for each class
                meanGauge{e} = nozzleGauge(filteredClassNozzle{e},{@(X,v)particleCounter(X,v),@(X,v)accumulatorFunc(X,@(V)sum(V,2),v)});
                % measure the average particle "pressure"
                M{e} = meanGauge{e}.takeMetric();
                % hook up covariance 
                covGauge{e} = nozzleGauge(filteredClassNozzle{e},{@(X,v)particleCOV(X,M{e},v)});
                % measure the covariance for the data stream
                C{e} = covGauge{e}.takeMetric();
                % store the whole data for each class if-and-only-if
                % there is a call for a mixture distribution
                if nGRPs(e) ~= 0
                    D{e} = filteredClassNozzle{e}.collectFullSpray();
                else 
                    D{e} = filteredClassNozzle{e}.collectFullSpray();
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % make hidden markov model
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [hmm] = makeChainRepeat(M{1}{2}/M{1}{1},C{1}{1},M{2}{2}/M{2}{1},C{2}{1},D,transPROB,repeatingUnits,nCOMP,nGRPs);




            % create transformation nozzle function that will apply the hidden
            % markov model
            func = @(X,e0,e1)(1-mod(hmm.Viterbi(X,ones(nCOMP,1),1),2));

        
            Y = [];
            
            
            %{
            cnt = 1;
            % reset the pointer to the reduced Image Nozzle
            reducedImageNozzle.resetPtr();
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get the data values
            while reducedImageNozzle.hasNext()
                % get the reduced image
                %tmpD{cnt} = reducedImageNozzle.next();
                cnt = cnt + 1;
                cnt
            end
            %}
            cnt = 1;
            %reducedImageNozzle.dataSource.dataSource.MEMstore = {};
            %for i = 1:(reducedImageNozzle.nozzlePtr.majorMax)
            while reducedImageNozzle.hasNext()
                %tmpD = reducedImageNozzle.read(i);
                tmpD = reducedImageNozzle.next();
                % apply the hidden-markov model to the signal
                tmpY = func(tmpD);
                trainD{cnt} = [tmpY;tmpD];
                % build up the new Y on this loop pass for re-cursive
                % training
                Y = [Y tmpY];
                % increment cnt
                cnt = cnt + 1;
                cnt
                %i
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% re-map the binary result from HMM application to the states of the 
            % HMM
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for c = 1:numel(trainD)
                % select the ones and the zeros
                sig1 = trainD{c}(1,:)==1;
                sig0 = trainD{c}(1,:)==0;
                % get the regions
                R1 = regionprops(sig1,'PixelIdxList');
                R0 = regionprops(sig0,'PixelIdxList');
                % make a zeros
                z = zeros(size(sig1));
                % make a state chain
                value1 = [2 4 6];
                value0 = [1 3 5 7];
                % for each region,set the new values for the HMM state chain
                for e = 1:numel(R1)
                    z(R1(e).PixelIdxList) = value1(e);
                end
                for e = 1:numel(R0)
                    z(R0(e).PixelIdxList) = value0(e);
                end
                % asssign the new values for the HMM state chain for
                % updating
                trainD{c}(1,:) = z;
            end

            %% call the update function
            hmm.update(trainD,ones(nCOMP,1));
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create transformation nozzle
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        func = @(X,e0,e1)(1-mod(hmm.Viterbi(X,ones(nCOMP,1),1),2));
        tester = dataNozzle(func,reducedImageNozzle,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    catch ME
        ME;
    end
    %{
    %%
    maskNozzle.resetPtr();
    tester.resetPtr();
    
    cnt = 1;
    while tester.hasNext()
        gidx = tester.next();
        tidx = maskNozzle.next();
        I = imageNozzle.read(cnt);
        I = reshape(I,[size(I,1)/3 3 size(I,2)]);
        I = permute(I,[1 3 2]);
        gidx = imresize(gidx,[1 size(I,2)]);
        tidx = imresize(tidx,[1 size(I,2)]);
        msk = repmat(gidx,[size(I,1) 1]);
        msk2 = repmat(tidx,[size(I,1) 1]);
        I = permute(I,[2 1 3]);
        out=  flattenMaskOverlay(double(I)/280,logical(msk'));
        out = flattenMaskOverlay(out,logical(msk2'),.1,'b');
        imshow(out,[]);
        drawnow
        %{
        plot(gidx,'r');
        hold on
        plot(tidx,'b')
        drawnow
        hold off
        %}
        cnt = cnt + 1;
    end
    %}
    %%
    %{
    tester.resetPtr();
    verticalANYnozzle.resetPtr();
    while tester.hasNext()
        gidx = tester.next();
        gidx = 1-mod(gidx,2);
        tidx = verticalANYnozzle.next();
        plot(gidx,'r')
        hold on
        plot(tidx,'b')
        drawnow
        hold off
    end
    %}
    
        %{
        % remove outiers for getting the number of groups
        perToR = .2;
        for e = 1:2
            toK{e} = 1:size(D{e},2);
            for d = 1:size(D{e},1)
                [fi,xi] = ksdensity(D{e}(d,:));
                fi = fi / sum(fi);
                fy = interp1(xi,fi,D{e}(d,:));
                [fys sidx] = sort(fy);
                toK{e} = intersect(toK{e},sidx(round(size(D{e},2)*perToR):end)');
            end
        end
        %}
        %{
        for e = 1:2
            if isnan(nGRPs(e))
                [IDX,CE,SUMD,K]=best_kmeans(D{e}(:,:)');
                if K > 8
                    K = K - 2;
                end
                nGRPs(e) = K;
            end
        end
        %}
    
        %{
        numComponents = [];
        AIC = [];
        RV = 0.00000000001;
        for e = 1:2
            for k = 1:10
                try
                    obj = fitgmdist(D{e}(:,toK{e})',k,'Replicates',3,'Start','plus','CovarianceType','diagonal','RegularizationValue',RV);
                    AIC(e,k)= obj.AIC + 2*k*(k+1)*(size(D{e},2)-k-1)^-1;
                catch
                    AIC(e,k) = inf;
                end
                plot(AIC(e,:))
                drawnow
            end
            [minAIC,numComponents(e)] = min(AIC(e,:));
            numComponents
        end
        %}
        
end