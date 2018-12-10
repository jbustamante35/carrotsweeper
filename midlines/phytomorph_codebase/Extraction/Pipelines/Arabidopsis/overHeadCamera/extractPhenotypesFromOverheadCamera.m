function [pHMeasure] = extractPhenotypesFromOverheadCamera(fileList,levelD,levelDD,oPath)
    try
        
        if isdeployed()
            %% remove dark images
            for e = 1:numel(fileList)
                I = imread(fileList{e});
                G = rgb2gray(I);
                value(e) = mean(G(:));
            end
            fileList(value < 50) = [];
        end
        
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
        
        
        
        mkdir(oPath)
        tMask = [];
        toOpSET = fliplr(fileList);
        size(toOpSET)
        nmVec = 1:numel(toOpSET);
        pHMeasure = [];
        close all
        for e = 1:numel(toOpSET)


            oI = double(imread(toOpSET{e}));
            osz = size(oI);
            I = permute(oI,[3 1 2]);
            sz = size(I);
            I = reshape(I,[sz(1) prod(sz(2:3))]);
            idx = cluster(levelD,double(I'));
            idx = reshape(idx,osz(1:2));
            %imshow(idx,[]);
            [newMask] = applyClusterStage(oI,idx==2,levelDD);
            plantMask = newMask == 2 | newMask == 1;
            plantMask0 = newMask == 2;
            plantMask1 = newMask == 1;

            plantMask = bwareaopen(plantMask,1000);
            %plantMask = imclearborder(plantMask);

            %{
            plantMask0 = bwareaopen(plantMask0,100);
            plantMask0 = imclearborder(plantMask0);

            plantMask1 = bwareaopen(plantMask1,100);
            plantMask1 = imclearborder(plantMask1);
            %}



            finalPlantMask = zeros(size(plantMask));

            %% blob removal

            if e == 1
                
                RRM = regionprops(logical(plantMask),'Eccentricity','PixelIdxList');
                plantMask = logical(zeros(size(plantMask)));
                for r = 1:numel(RRM)
                    if RRM(r).Eccentricity < .9
                        plantMask(RRM(r).PixelIdxList) = 1;
                    end
                end
                
                
            end


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
            out = flattenMaskOverlay(out,plantMask1,.5,'b');


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



            saveas(gca,[oPath num2str(nmVec(e)) '_overlay.jpg'])




        end
        pHMeasure.Area = fliplr(pHMeasure.Area);
        pHMeasure.Area0 = fliplr(pHMeasure.Area0);
        pHMeasure.Area1 = fliplr(pHMeasure.Area1);
        pHMeasure.ConvexArea = fliplr(pHMeasure.ConvexArea);
        pHMeasure.MajorAxisLength = fliplr(pHMeasure.MajorAxisLength);
        pHMeasure.MinorAxisLength = fliplr(pHMeasure.MinorAxisLength);

        csvwrite([oPath 'phenotype_Area.csv'],pHMeasure.Area)
        csvwrite([oPath 'phenotype_Area0.csv'],pHMeasure.Area0)
        csvwrite([oPath 'phenotype_Area1.csv'],pHMeasure.Area1)
        csvwrite([oPath 'phenotype_ConvexArea.csv'],pHMeasure.ConvexArea)
        csvwrite([oPath 'phenotype_MajorAxisLength.csv'],pHMeasure.MajorAxisLength)
        csvwrite([oPath 'phenotype_MinorAxisLength.csv'],pHMeasure.MinorAxisLength)
        close all
    catch ME
        getReport(ME)
        close all
    end
end