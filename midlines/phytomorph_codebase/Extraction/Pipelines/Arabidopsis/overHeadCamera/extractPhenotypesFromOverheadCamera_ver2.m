function [pHMeasure] = extractPhenotypesFromOverheadCamera_ver2(fileList,cT,oPath)
    try
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isdeployed()
            %% remove dark images
            for e = 1:numel(fileList)
                I = imread(fileList{e});
                G = rgb2gray(I);
                value(e) = mean(G(:));
            end
            fileList(value < 50) = [];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sort data
        for e = 1:numel(fileList)
            [tp,tn,te] = fileparts(fileList{e});
            newN = [];
            for n = 1:numel(tn)
                nn = str2num(tn(n));
                if ~isempty(nn)
                    newN = [newN tn(n)];
                end
            end
            %fidx = strfind(tn,'_');
            %nm(e)=  str2num(tn((fidx(1)+1):end));
            nm(e) = str2num(newN);
        end
        [~,sidx] = sort(nm);
        fileList = fileList(sidx);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmpU = imread(fileList{1});
        tmpU = zeros(size(tmpU));
        
        for e = 1:numel(fileList)
            tmpU = tmpU + double(imread(fileList{e}))/255;
        end
        tmpU = tmpU *numel(fileList)^-1;
        initMask = cT.clusterImage(tmpU*255);
        initMask = (initMask == 2);
        initMask = logical(bwareaopen(initMask,400));
        
        
        clusterR = regionprops(logical(initMask),'Centroid','PixelIdxList');
        centroids = [];
        for r = 1:numel(clusterR)
            centroids = [centroids;clusterR(r).Centroid];
        end
        Z = linkage(centroids,'centroid','euclidean');
        plantBlobClusters = cluster(Z,'cutoff',200,'criterion','distance');
        
        
        UQ = unique(plantBlobClusters);
        plantCenterPoint = [];
        plantLabelMatrix = zeros(size(initMask));
        for b = 1:numel(UQ)
            fidx = find(plantBlobClusters == UQ(b));
            tmpClusterCenters = [];
            for p = 1:numel(fidx)
                tmpClusterCenters = [tmpClusterCenters;clusterR(fidx(p)).Centroid];
                plantLabelMatrix(clusterR(fidx(p)).PixelIdxList) = b;
            end
            plantCenterPoint = [plantCenterPoint ; mean(tmpClusterCenters,1)];
        end
        %{
        R = regionprops(initMask,'EquivDiameter','Centroid');
        plantCenterPoint = [];
        for r = 1:numel(R)
            plantCenterPoint = [plantCenterPoint;R(r).Centroid];
        end
        %}
        
        
        DM = squareform(pdist(plantCenterPoint));
        for d = 1:size(DM,1)
            DM(d,d) = inf;
        end
        RAD = mean(min(DM,[],1));
        RAD = RAD*.7;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        mkdir(oPath)
        tMask = [];
        toOpSET = fliplr(fileList);
        size(toOpSET)
        nmVec = 1:numel(toOpSET);
        pHMeasure = [];
        close all
        
        
        for e = 1:numel(toOpSET)

            newMask = cT.clusterImage(toOpSET{e});
            oI = double(imread(toOpSET{e}));
            osz = size(oI);
            %{
            I = permute(oI,[3 1 2]);
            sz = size(I);
            I = reshape(I,[sz(1) prod(sz(2:3))]);
            
            
            idx = cluster(levelD,double(I'));
            idx = reshape(idx,osz(1:2));
            
            
            
            %imshow(idx,[]);
            %}
            %[newMask] = applyClusterStage(oI,idx==2,levelDD);
            
            %plantMask = newMask == 2 | newMask == 1;
            %plantMask0 = newMask == 2;
            %plantMask1 = newMask == 1;
            
            
            
            plantMask = newMask == 2;
           
            %plantMask0 = plantMask;
            %plantMask1 = plantMask;
            %plantMask = imclearborder(plantMask);

            %{
            plantMask0 = bwareaopen(plantMask0,100);
            plantMask0 = imclearborder(plantMask0);

            plantMask1 = bwareaopen(plantMask1,100);
            plantMask1 = imclearborder(plantMask1);
            %}



            plantMask = bwareaopen(logical(plantMask),400);

            %% blob removal
            %{
            if e == 1
                plantMask = bwareaopen(logical(plantMask),1000);
                % remove objects based on length width ratio
                RRM = regionprops(logical(plantMask),'Eccentricity','PixelIdxList');
                plantMask = logical(zeros(size(plantMask)));
                for r = 1:numel(RRM)
                    if RRM(r).Eccentricity < .90
                        plantMask(RRM(r).PixelIdxList) = 1;
                    end
                end
                
                
                clusterR = regionprops(logical(plantMask),'Centroid','PixelIdxList');
                centroids = [];
                for r = 1:numel(clusterR)
                    centroids = [centroids;clusterR(r).Centroid];
                end
                Z = linkage(centroids,'centroid','euclidean');
                plantBlobClusters = cluster(Z,'cutoff',200,'criterion','distance');


                UQ = unique(plantBlobClusters);
                plantCenterPoint = [];
                plantLabelMatrix = zeros(size(plantMask));
                for b = 1:numel(UQ)
                    fidx = find(plantBlobClusters == UQ(b));
                    tmpClusterCenters = [];
                    for p = 1:numel(fidx)
                        tmpClusterCenters = [tmpClusterCenters;clusterR(fidx(p)).Centroid];
                        plantLabelMatrix(clusterR(fidx(p)).PixelIdxList) = b;
                    end
                    plantCenterPoint = [plantCenterPoint ; mean(tmpClusterCenters,1)];
                end
                RGB_label = label2rgb(plantLabelMatrix);
            end
            %}

            
            finalPlantMask = zeros(size(plantMask));
            RPM = regionprops(logical(plantMask),'Eccentricity','PixelIdxList','Area','Centroid');
            
            % for 
            for p = 1:size(plantCenterPoint,1)
                tmpD = [];
                for b = 1:numel(RPM)
                    % get working plant centroid
                    tmpCen = RPM(b).Centroid;
                    % if with radius
                    tmpD(b) = norm(tmpCen - plantCenterPoint(p,:));
                    if tmpD(b) < RAD
                        finalPlantMask(RPM(b).PixelIdxList) = 1;
                    end
                end
                fidx = find(tmpD < RAD);
                pHMeasure.Area(p,e) = sum([RPM(fidx).Area]);
            end
            
            
            out = flattenMaskOverlay(oI/255,logical(finalPlantMask),.5,'r');
            close all
            image(out);
            viscircles(plantCenterPoint,RAD*ones(size(plantCenterPoint,1),1));
            axis off
            axis equal
            drawnow
            saveas(gca,[oPath num2str(nmVec(e)) '_overlay.jpg']);
            hold on
            for p = 1:size(plantCenterPoint,1)
                text(plantCenterPoint(p,1),plantCenterPoint(p,2),num2str(p));
            end
            saveas(gca,[oPath num2str(nmVec(e)) '_overlayLabels.jpg']);
            
            
            
            %{
            
            CP = [];
            BOX = {};
            CONN = {};
            if e == 1
                pHMeasure.Area = [];
                Rinit = regionprops(logical(plantMask),'PixelIdxList','Area','Centroid','ConvexHull','BoundingBox','EquivDiameter','MajorAxisLength','MinorAxisLength','ConvexArea');
                fprintf(['Found ' num2str(numel(Rinit)) ' plants in last image \n']);
                for r = 1:numel(Rinit)
                    pHMeasure.Area(r,e) = Rinit(r).Area;
                    %text(Rinit(r).Centroid(1),Rinit(r).Centroid(2),num2str(r));
                    CP(r,:) = Rinit(r).Centroid;
                    CONN{r} = Rinit(r).ConvexHull;
                    BOX{r} = Rinit(r).BoundingBox;
                    pHMeasure.ConvexArea(r,e) = Rinit(r).ConvexArea;
                    pHMeasure.MajorAxisLength(r,e) = Rinit(r).MajorAxisLength;
                    pHMeasure.MinorAxisLength(r,e) = Rinit(r).MinorAxisLength;
                    finalPlantMask(Rinit(r).PixelIdxList) = 1;
                end
            else
                R = regionprops(logical(plantMask),'PixelIdxList','Area','Centroid','ConvexHull','BoundingBox','EquivDiameter','MajorAxisLength','MinorAxisLength','ConvexArea');
                fprintf(['Found ' num2str(numel(R)) ' plants in ' num2str(e) ' image \n']);
                   
                for ri = 1:numel(Rinit)
                    tmpM = zeros(size(plantMask));

                    dist = [];
                    for r = 1:numel(R)
                        dist(r) = norm(Rinit(ri).Centroid - R(r).Centroid);
                        %{
                        if ~isempty(intersect(Rinit(ri).PixelIdxList,R(r).PixelIdxList))
                            tmpM(R(r).PixelIdxList) = 1;
                        end
                        %}
                    end
                    [dist,midx] = min(dist);
                    tmpM(R(midx).PixelIdxList) = 1;

                    
                    Rtmp = regionprops(logical(tmpM),'PixelIdxList','Area','Centroid','ConvexHull','BoundingBox','EquivDiameter','MajorAxisLength','MinorAxisLength','ConvexArea');
                    fprintf(['Found ' num2str(numel(Rtmp)) ' plants in ' num2str(e) ' image \n']);
                    
                    
                    % select first for now:select closes later

                    pHMeasure.Area(ri,e) = Rtmp(1).Area;
                    %text(R(r).Centroid(1),R(r).Centroid(2),num2str(ri));
                    CP(ri,:) = Rtmp(1).Centroid;
                    CONN{ri} = Rtmp(1).ConvexHull;
                    BOX{ri} = Rtmp(1).BoundingBox;
                    pHMeasure.ConvexArea(ri,e) = Rtmp(1).ConvexArea;
                    pHMeasure.MajorAxisLength(ri,e) = Rtmp(1).MajorAxisLength;
                    pHMeasure.MinorAxisLength(ri,e) = Rtmp(1).MinorAxisLength;
                    finalPlantMask(Rtmp(1).PixelIdxList) = 1;
                end
            end

            plantMask0 = logical(plantMask0.*finalPlantMask);
            plantMask1 = logical(plantMask1.*finalPlantMask);

            %tMask(:,:,:,e) = cat(3,plantMask0,plantMask1,plantMask,finalPlantMask);


            out = flattenMaskOverlay(oI/255,plantMask0,.5,'r');
            %out = flattenMaskOverlay(out,plantMask1,.5,'b');


            close all
            image(out)
            hold on

            for r = 1:size(CP,1)
                text(CP(r,1)-50,CP(r,2)-50,num2str(r),'Color','y');
                plot(CONN{r}(:,1),CONN{r}(:,2),'r');
                tmpI0 = imcrop(plantMask0,BOX{r});
                tmpI1 = imcrop(plantMask1,BOX{r});

                pHMeasure.Area0(r,e) = sum(tmpI0(:));
                pHMeasure.Area1(r,e) = sum(tmpI1(:));




            end
            axis off

            %}

           

            close all

        end
        pHMeasure.Area = fliplr(pHMeasure.Area);
        %{
        pHMeasure.Area0 = fliplr(pHMeasure.Area0);
        pHMeasure.Area1 = fliplr(pHMeasure.Area1);
        pHMeasure.ConvexArea = fliplr(pHMeasure.ConvexArea);
        pHMeasure.MajorAxisLength = fliplr(pHMeasure.MajorAxisLength);
        pHMeasure.MinorAxisLength = fliplr(pHMeasure.MinorAxisLength);
        %}
       
        csvwrite([oPath 'phenotype_Area.csv'],pHMeasure.Area)
        %{
        csvwrite([oPath 'phenotype_Area0.csv'],pHMeasure.Area0)
        csvwrite([oPath 'phenotype_Area1.csv'],pHMeasure.Area1)
        csvwrite([oPath 'phenotype_ConvexArea.csv'],pHMeasure.ConvexArea)
        csvwrite([oPath 'phenotype_MajorAxisLength.csv'],pHMeasure.MajorAxisLength)
        csvwrite([oPath 'phenotype_MinorAxisLength.csv'],pHMeasure.MinorAxisLength)
        %}
        close all
    catch ME
        getReport(ME)
        close all
    end
end