function [grade,T,BO,out] = opti3(source,target,para,oI,initDelta,paraLabels)
    try
        T = 0;
        if numel(target) == 1
            T = target;
            target = [];
        end
        grade = NaN;

        init = initDelta(paraLabels==2);
        delta = initDelta(paraLabels==3);

        para = init + (para - delta);
        % para(1) = marker subtraction
        % para(2:4) = wieghts for mean

        % para(5) = lower limit for objects for inner
        % para(6) = upper limit for objects for inner
        % para(7) = abs prob threshold for inner
        % para(8) = relative threshold for inner

        % para(9) = lower limit for objects for inner
        % para(10) = upper limit for objects for inner
        % para(11) = abs prob threshold for inner
        % para(12) = relative threshold for inner

        % para(13) = smooth
        ag = [];
        ret1 = [];
        BO = [];
        disp = false;

        parfor o = 1:size(source,4)
            % assign the source to tmp
            tmpS = source(:,:,:,o);
            % perform custom weighting
            for slice = 1:size(tmpS,3)
                tmpS(:,:,slice) = tmpS(:,:,slice) * para(slice+1);
            end
            % perform filter smoothing
            tmpS = imfilter(tmpS,fspecial('gaussian',[21 21],max(round(para(13)),1)),'replicate');
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


            POT = double(BO).*wSUM;

            BBO = imclearborder(BO) == 0 & BO == 1;
            BO = imclearborder(BO);

            R_inner_objects_abs = regionprops(BO,POT,'MeanIntensity','Area','PixelIdxList','Centroid','MajorAxisLength');
            R_inner_objects_rel = regionprops(BO,PEAKS,'MeanIntensity','Area','PixelIdxList','Centroid','MajorAxisLength');
            R_edge_objects_abs = regionprops(BBO,POT,'MeanIntensity','Area','PixelIdxList');
            R_edge_objects_rel = regionprops(BBO,PEAKS,'MeanIntensity','Area','PixelIdxList');

            if ~isempty(target)
                tmpT = logical(target(:,:,o));
                tmpT = imerode(tmpT,strel('disk',2,0));
                Rtmp = regionprops(tmpT,'PixelIdxList');
                VAL{o} = zeros(numel(R_inner_objects_abs),1);
                for p = 1:numel(R_inner_objects_abs)
                    for q = 1:numel(Rtmp)
                        numI = intersect(Rtmp(q).PixelIdxList,R_inner_objects_abs(p).PixelIdxList);
                        if numI > 20
                            VAL{o}(p) = 1;
                        end
                    end
                end
            end


            %{
            Features_inner = [[R_inner_objects_abs.Area]' [R_inner_objects_abs.MeanIntensity]' ...
                        [R_inner_objects_rel.Area]' [R_inner_objects_rel.MeanIntensity]'];
            Features_edge = [[R_edge_objects_abs.Area]' [R_edge_objects_abs.MeanIntensity]' ...
                        [R_edge_objects_rel.Area]' [R_edge_objects_rel.MeanIntensity]'];
            %}
            Features_inner = [[R_inner_objects_abs.Area]' [R_inner_objects_abs.MeanIntensity]' ...
                        [R_inner_objects_rel.Area]' ([R_inner_objects_abs.Area].*[R_inner_objects_abs.MeanIntensity])'];
            Features_edge = [[R_edge_objects_abs.Area]' [R_edge_objects_abs.MeanIntensity]' ...
                        [R_edge_objects_rel.Area]' ([R_edge_objects_abs.Area].*[R_edge_objects_abs.MeanIntensity])'];
            Features_extra_inner = [[R_inner_objects_abs.MajorAxisLength]'];
                    
            if ~isempty(Features_inner)
                REPO{o} = [Features_inner Features_extra_inner];
            else
                REPO{o} = single(NaN*ones(1,5));
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
                    
            if ~isempty(Features_inner)
                innerIDX = Features_inner(:,1) > para(5) & ...
                Features_inner(:,1) < para(6) & ...
                Features_inner(:,2) > para(7) & ...
                Features_inner(:,4) > para(8);
                filter{o} = double(innerIDX);
                R_inner_objects_abs = R_inner_objects_abs(innerIDX);
            else
                filter{o} = NaN;
            end
            
            if ~isempty(Features_edge)
                edgeIDX = Features_edge(:,1) > para(9) & ...
                Features_edge(:,1) < para(10) & ...
                Features_edge(:,2) > para(11) & ...
                Features_edge(:,4) > para(12);
                filter_E{o} = double(edgeIDX);
                R_edge_objects_abs = R_edge_objects_abs(edgeIDX);
            else
                filter_E{o} = NaN;
            end


            
            %{  
                figure
                imshow(out,[])
                hold on;
                for p = 1:numel(R_inner_objects_abs)
                    text(R_inner_objects_abs(p).Centroid(1),R_inner_objects_abs(p).Centroid(2),num2str(p),'Color','r');
                end
            %}


            % repaint
            Z = zeros(size(oI,1),size(oI,2));

            for q1 = 1:numel(R_inner_objects_abs)
                Z(R_inner_objects_abs(q1).PixelIdxList) = 1;
            end
            for q1 = 1:numel(R_edge_objects_abs)
                Z(R_edge_objects_abs(q1).PixelIdxList) = 1;
            end

            outCount(o) = numel(R_inner_objects_abs) + numel(R_edge_objects_abs);


            hS{o} = Z;

            if ~isempty(target)

                tmpT = logical(target(:,:,o));
                tmpT = imerode(tmpT,strel('disk',2,0));
                Rtmp = regionprops(tmpT);

                targetCount(o) = numel(Rtmp);

                targets = [tmpT(:) ~tmpT(:)];
                outputs = [Z(:) ~Z(:)];


                TOTX{o} = targets(:,1);
                TOTY{o} = outputs(:,1);
            end


        end

        if ~isempty(target)

            


            STORE = TOTY;

            VAL = cell2mat(VAL');
            REPO = cell2mat(REPO');
            REPO_E = cell2mat(REPO_E');
            REPO(find(any(isnan(REPO),2)),:) = [];
            REPO_E(find(any(isnan(REPO_E),2)),:) = [];


            Y = ind2vec(VAL'+1);
            X = REPO';
            PT = patternnet(5);
            PT = train(PT,X,Y);

            %{
            filter = find(REPO(:,1) > para(5) & REPO(:,1) < para(6));
            filter_E = find(REPO_E(:,1) > para(9) & REPO_E(:,1) < para(10));
            %}
            filter = cell2mat(filter');
            filter_E = cell2mat(filter_E');
            filter(isnan(filter)) = [];
            filter = logical(filter);
            filter_E(isnan(filter_E)) = [];
            filter_E = logical(filter_E);


            cidx1 = count(REPO(filter,1));
            [cidx2,value] = count(REPO(filter,5));
            cidx1_E = count(REPO_E(filter_E,1));



            nC2 = zeros(size(REPO,1),1);
            nC = zeros(size(REPO,1),1);
            if sum(filter) > 0
                nC(filter) = cidx1;
                nC2(filter) = cidx2;
            end
            nC_E = zeros(size(REPO_E,1),1);
            if sum(filter_E) > 0
                nC_E(filter_E) = cidx1_E;
            end
            


            SUBS = cell2mat(SUBS');
            SUBS_E = cell2mat(SUBS_E');

            fidx = find(cidx1~=0);
            altStack = [];

           


            for ni = 1:size(source,4)
                tmpZ = zeros(size(oI,1),size(oI,2));
                tmpZ2 = zeros(size(oI,1),size(oI,2));
                sidx = find(SUBS(:,1) == ni & (nC==1));
                for s = 1:numel(sidx)
                    IDX = SUBS(sidx(s),:);
                    tmpZ(REPO_INNER_ABS{IDX(1)}(IDX(2)).PixelIdxList) = 1;
                end

                sidx2 = find(SUBS_E(:,1) == ni & (nC_E==1));
                for s = 1:numel(sidx2)
                    IDX = SUBS_E(sidx2(s),:);
                    tmpZ(REPO_EDGE_ABS{IDX(1)}(IDX(2)).PixelIdxList) = 1;
                end

                sidx3 = find(SUBS(:,1) == ni & nC2==2);
                for s = 1:numel(sidx3)
                    IDX = SUBS(sidx3(s),:);
                    tmpZ2(REPO_INNER_ABS{IDX(1)}(IDX(2)).PixelIdxList) = 1;
                    %imshow(tmpZ2,[]);
                    %drawnow
                end

                altSTORE{ni} = tmpZ;
                altSTORE2{ni} = tmpZ2;
                altStack = [altStack;tmpZ(:)+tmpZ2(:)];
                twoCount(ni) = numel(sidx3);
                newCount(ni) = numel(sidx) + numel(sidx2) + numel(sidx3);
            end
            
            
            
            TOTX = cell2mat(TOTX');
            TOTY = cell2mat(TOTY');
            %{
            [Xp,Yp,T,grade,OPTROCPT] = perfcurve(TOTX,TOTY,1);
            T = T(Xp==OPTROCPT(1) & Yp==OPTROCPT(2));
            T = 0;
            %}
            T = 0;
            %[tpr,fpr,thresholds] = roc(logical(TOTX'),[logical(TOTY');~(TOTY'>0)]);
            [grade] = matthews_correlation(logical(TOTX),logical(TOTY));
            [grade2] = matthews_correlation(logical(TOTX),logical(altStack));
            %grade = tpr(2)*fpr(2)^-1;
            SSS = {83 19 76 75 34 15 22 23 98 30 31 32 42 43 44 51 52 53};

            MASTER = [];
            MASTERC = [];
            LAB = '';
            for wow = 1:numel(SSS)
                sel = SSS{wow};
                toShow = reshape(STORE{sel},[size(target,1),size(target,2)]);
                toShow2 = reshape(altSTORE{sel},[size(target,1),size(target,2)]);
                toShow3 = reshape(altSTORE2{sel},[size(target,1),size(target,2)]);
                out = flattenMaskOverlay(oI(:,:,sel),toShow > T,.5,'r');
                out2 = flattenMaskOverlay(oI(:,:,sel),toShow2 > T,.5,'r');
                out2 = flattenMaskOverlay(out2,toShow3 > T,.5,'g');
                if ~isempty(target)
                    out = flattenMaskOverlay(out,logical(target(:,:,sel)),.5,'b');
                    out2 = flattenMaskOverlay(out2,logical(target(:,:,sel)),.5,'b');
                end
                
               MASTERC = cat(1,MASTERC,out2);
               if mod(wow,3) == 0
                MASTER = [MASTER MASTERC];
                MASTERC = [];
                end
                if wow == 1
                    LAB = [num2str(newCount(sel)) '-->' num2str(targetCount(sel))];
                end
            end


            imshow(MASTER,[]);
            title(LAB);
            drawnow




            grade = -grade;
            grade2 = -grade2;

            grade(2) = -corr(targetCount',outCount');
            grade2(2) = -corr(targetCount',newCount');
            grade(3) = .25*mean(abs(outCount - targetCount));
            grade2(3) = .25*mean(abs(newCount - targetCount));
            if sum(newCount) == 329 & sum(outCount) == 329
            here = 1;
            end
            [sum(newCount) sum(outCount) sum(targetCount)];
           
            grade(isnan(grade)) = 100;
            grade2(isnan(grade2)) = 100;
            grade;
            grade2
            grade = sum(grade);
            grade = sum(grade2);

            if any(twoCount==2)
                tidx = find(twoCount==2);
                if any(targetCount(tidx) == newCount(tidx))
                    stop = 1;
                end
            end

            %{
            imshow(ret2(:,:,:,1),[]);
            drawnow
            %waitforbuttonpress
            %-(mean(grade))
            %grade = -(mean(grade));
            % grade;
            rmidx = isinf(gr2);
            gr2(rmidx) = NaN;
            MX = max(gr2);
            gr2(isnan(gr2)) = MX;
            %grade = -mean(gr);
            grade = -mean(gr2);
            if grade < -100
                here = 1;
            end
            %}
        else
            BO = hS{1} > 0;
            out = flattenMaskOverlay(oI(:,:,1),BO,.5,'r');
        end


        para;
    catch ME
        ME
    end
    
    
    
    
        
        
end