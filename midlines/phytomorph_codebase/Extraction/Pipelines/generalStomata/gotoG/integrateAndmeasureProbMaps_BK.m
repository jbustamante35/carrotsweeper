function [grade,PT,PT2,mask,out] = integrateAndmeasureProbMaps(source,target,para,oI,initDelta,paraLabels,PT,PT2)
    try
        % para(1:3) := smooth parameter
        % para(4) : = smooth parameter
        % para(5) : = reconstruction para
        T = 0;
        mask = [];
        out = [];
        grade = NaN;
        REPO_INNER_ABS = {};
        REPO_EDGE_ABS = {};
        REPO = {};
        REPO_E = {};
        SUBS = {};
        VAL = {};



        init = initDelta(paraLabels==2);
        delta = initDelta(paraLabels==3);
        para = init + (para - delta);

        for o = 1:size(source,4)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % assign the source to tmp
            tmpS = source(:,:,:,o);
            % perform custom weighting
            for slice = 1:size(tmpS,3)
                tmpS(:,:,slice) = tmpS(:,:,slice) * para(slice+1);
            end
            % perform filter smoothing
            tmpS = imfilter(tmpS,fspecial('gaussian',[21 21],max(round(para(5)),1)),'replicate');
            % take the mean of weighted images
            wSUM = mean(tmpS,3);
            % perform subtract for marker and mask
            LEVEL = para(1);
            sumMarker = wSUM - LEVEL;
            RECON = imreconstruct(sumMarker,wSUM);
            PEAKS = wSUM - RECON;
            BO = PEAKS > 0;
            % find the location which is above the reconstruction
            BO = imfill(BO,'holes');
            BO = bwareaopen(BO,20);
            % fuse together
            BO = imclose(BO,strel('disk',3,0));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            POT = double(BO).*wSUM;
            BBO = imclearborder(BO) == 0 & BO == 1;
            %BO = imclearborder(BO);

            R_inner_objects_abs = regionprops(BO,POT,'MeanIntensity','Area','PixelIdxList','Centroid','MajorAxisLength');
            R_inner_objects_rel = regionprops(BO,PEAKS,'MeanIntensity','Area','PixelIdxList','Centroid','MajorAxisLength');
            R_edge_objects_abs = regionprops(BBO,POT,'MeanIntensity','Area','PixelIdxList');
            R_edge_objects_rel = regionprops(BBO,PEAKS,'MeanIntensity','Area','PixelIdxList');

            %{
            Features_inner = [[R_inner_objects_abs.Area]' [R_inner_objects_abs.MeanIntensity]' ...
                        [R_inner_objects_rel.Area]' [R_inner_objects_rel.MeanIntensity]'];
            Features_edge = [[R_edge_objects_abs.Area]' [R_edge_objects_abs.MeanIntensity]' ...
                        [R_edge_objects_rel.Area]' [R_edge_objects_rel.MeanIntensity]'];
            %}
            Features_inner = [[R_inner_objects_abs.Area]' [R_inner_objects_abs.MeanIntensity]' ...
                             ([R_inner_objects_abs.Area].*[R_inner_objects_abs.MeanIntensity])' ...
                              [R_inner_objects_rel.MeanIntensity]' ...
                             ([R_inner_objects_rel.Area].*[R_inner_objects_rel.MeanIntensity])'];
            Features_edge = [[R_edge_objects_abs.Area]' [R_edge_objects_abs.MeanIntensity]' ...
                        [R_edge_objects_rel.Area]' ([R_edge_objects_abs.Area].*[R_edge_objects_abs.MeanIntensity])'];
            Features_extra_inner = [[R_inner_objects_abs.MajorAxisLength]'];


            if ~isempty(Features_inner)
                REPO{o} = [Features_inner Features_extra_inner];
            else
                REPO{o} = single(NaN*ones(1,6));
            end
            if ~isempty(Features_edge)
                REPO_E{o} = Features_edge;
            else
                REPO_E{o} = single(NaN*ones(1,4));
            end




            SUBS{o} = [o*ones(numel(R_inner_objects_abs),1) [1:numel(R_inner_objects_abs)]'];
            SUBS_E{o} = [o*ones(numel(R_edge_objects_abs),1) [1:numel(R_edge_objects_abs)]'];
            REPO_INNER_ABS{o} = R_inner_objects_abs;
            REPO_EDGE_ABS{o} = R_edge_objects_abs;
            spotSTACK(:,:,o) = BO;



            if ~isempty(target)
                overlapThreshold = 20;
                tmpT = logical(target(:,:,o));
                tmpT = imerode(tmpT,strel('disk',2,0));
                Rtmp = regionprops(tmpT,'PixelIdxList');
                VAL{o} = zeros(numel(R_inner_objects_abs),1);
                targetCount(o) = numel(Rtmp);
                for p = 1:numel(R_inner_objects_abs)
                    for q = 1:numel(Rtmp)
                        numI = intersect(Rtmp(q).PixelIdxList,R_inner_objects_abs(p).PixelIdxList);
                        if numI > overlapThreshold
                            VAL{o}(p) = 1;
                        end
                    end
                end
            end



        end






        if ~isempty(target)


            VAL = cell2mat(VAL');
            REPO = cell2mat(REPO');
            rmidx = find(any(isnan(REPO),2));
            REPO(rmidx,:) = [];
            SUBS = cell2mat(SUBS');
            if ~isempty(VAL)

                PT = patternnet(10);
                Y = ind2vec(VAL'+1,2);
                X = REPO';
                PT.trainParam.showCommandLine = false;
                PT.trainParam.showWindow = false;
                PT = train(PT,X,Y);
                
                Ypre = PT(X);
                Ypre = (Ypre(2,:) > .5)';
            else
                Ypre = [];
            end



            ZSTACK = [];
            ASTACK = [];
            for t = 1:size(source,4)
                tmpZ = zeros(size(oI,1),size(oI,2));
                altZ = zeros(size(oI,1)*size(oI,2),size(REPO,2));
                sidx = find(SUBS(:,1)==t);

                for s = 1:numel(sidx)
                    if Ypre(sidx(s)) == 1
                        tmpZ(REPO_INNER_ABS{SUBS(sidx(s),1)}(SUBS(sidx(s),2)).PixelIdxList) = 1;
                    end
                    nn = size(REPO_INNER_ABS{SUBS(sidx(s),1)}(SUBS(sidx(s),2)).PixelIdxList,1);
                    altZ(REPO_INNER_ABS{SUBS(sidx(s),1)}(SUBS(sidx(s),2)).PixelIdxList,:) = repmat(REPO(sidx(s),:),[nn 1]);
                end
                REG = regionprops(logical(tmpZ));
                %altZ = reshape(altZ,[108 108 size(REPO,2)]);

                outCount(t) = numel(REG);
                ZSTACK(:,:,t) = tmpZ;
                ASTACK = [ASTACK;altZ];
            end



            TOTY = ZSTACK(:);
            TOTX = target(:);
            
            fidx1 = find(TOTX==1);
            fidx0 = find(TOTX==0);

            fidx0 = fidx0(randperm(numel(fidx0)));
            fidxM = [fidx1;fidx0(1:4*numel(fidx1))];

            PT2 = patternnet(10);
            Y = full(ind2vec(TOTX'+1,2));
            X = ASTACK';
            PT2.trainParam.showCommandLine = false;
            PT2.trainParam.showWindow = false;
            PT2 = train(PT2,X(:,fidxM),Y(:,fidxM),'useParallel','yes');
            Ypre2 = PT2(X);
            Ypre2 = (Ypre2(2,:) > .5)';
            

            YpreS = reshape(Ypre2,[108 108 100]);
            for oh = 1:size(YpreS,3)
                R = regionprops(logical(YpreS(:,:,oh)));
                newCount(oh) = numel(R);
            end

            TOTY = Ypre2;
            [grade] = -matthews_correlation(logical(TOTX),logical(TOTY));
            grade(2) = -corr(targetCount',newCount');
            grade(3) = .25*mean(abs(outCount - targetCount));
            grade(isnan(grade)) = 100;
            grade = sum(grade);
            grade

        end


            %{
            SSS = {83 19 76 75 34 15 22 23 98 30 31 32 42 43 44 51 52 53};
            MASTER = [];
            MASTERC = [];
            LAB = '';
            for wow = 1:numel(SSS)
                sel = SSS{wow};
                toShow = ZSTACK(:,:,sel);
                out = flattenMaskOverlay(oI(:,:,sel),toShow > T,.5,'r');
                if ~isempty(target)
                    out = flattenMaskOverlay(out,logical(target(:,:,sel)),.5,'b');
                end

               MASTERC = cat(1,MASTERC,out);
               if mod(wow,3) == 0
                MASTER = [MASTER MASTERC];
                MASTERC = [];
                end
                if wow == 1
                    LAB = [num2str(outCount(sel)) '-->' num2str(targetCount(sel))];
                end
            end
            imshow(MASTER,[]);
            title(LAB);
            drawnow
            %}

    catch ME
        ME
    end






end