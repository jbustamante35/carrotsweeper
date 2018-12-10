function [] = singleSeedlingImage(imageFile,smoothValue,threshSIG,EXT,topTRIM,SNIP,BKBOUND,oPath,rPath)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % imageFile : the image file to operate on
    % smoothValue : the smooth value for the integration signal along 1 dim
    % threshSIG : the threshold for finding the cone-tainers
    % EXT : the extension around the cone-tainers left and right
    % topTRIM : the amount to trim off the top
    % SNIP : the amont to use to find the base of the plant    
    % oPath : the location to save the results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % remove double quotes
    %{
    imageFile = removeDoubleQuotes(imageFile);
    smoothValue = removeDoubleQuotes(smoothValue);
    threshSIG = removeDoubleQuotes(threshSIG);
    EXT = removeDoubleQuotes(EXT);
    topTRIM = removeDoubleQuotes(topTRIM);
    SNIP = removeDoubleQuotes(SNIP);
    oPath = removeDoubleQuotes(oPath);
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert the strings to numbers if they are strings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ischar(smoothValue)
        smoothValue = str2num(smoothValue);
    end
    if ischar(threshSIG)
        threshSIG = str2num(threshSIG);
    end
    if ischar(EXT)
        EXT = str2num(EXT);
    end
    if ischar(topTRIM)
        topTRIM = str2num(topTRIM);
    end
    if ischar(SNIP)
        SNIP = str2num(SNIP);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % init the compute environment
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % init the icommands and create output directory
    %initIrods();
    mkdir(oPath);
    pushList = {};
    if isdeployed
        javaaddpath([pwd filesep 'core-3.2.1.jar']);
        javaaddpath([pwd filesep 'javase-3.2.1.jar']);
    end
    try
        OFFSET = 20;
        sigFILL = 1100;
        eT = 120;
        thresP = .1;
        TOP_THRESH = 1150;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % load the image, make gray, edge 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['starting: image load, gray, and edge \n']);
        fprintf(['working on: ' imageFile '\n']);
        % path and name
        [p nm] = fileparts(imageFile);
        % read the image
        I = imread(imageFile);
        % rectifiy
        [I angle] = rectifyImage(I);
        % get QR code
        nm = getQRcode(I);
        % crop off qr code
        I(1:TOP_THRESH,:,:) = [];
        % trim off X pixles from rotation
        I(:,1:30,:) = [];
        I(:,end-70:end,:) = [];
        % make gray scale
        G = rgb2gray(I);
        % filter the image
        G = imfilter(G,fspecial('gaussian',[13 13],4),'replicate');
        % find edge
        E = edge(G);
        fprintf(['ending: image load, gray, and edge \n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % load the image, make gray, edge 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find the cone-tainers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['starting: find the cone-tainer \n']);
        G = rgb2gray(single(I)/255);
        % smooth the image
        G = imfilter(G,fspecial('gaussian',[11 11]));
        % take grad
        [d1 d2] = gradient(G);
        % threshold grad along vertical direction
        d2 = abs(d2) > graythresh(abs(d2));
        %
        d2 = bwareaopen(d2,50);
        d2 = imclose(d2,strel('disk',10));
        sig = sum(abs(d2),1);
        sig = imfilter(sig,fspecial('average',[1 smoothValue]),'replicate');
        sig = bindVec(sig);
        threshSIG = graythresh(sig);
        % find the gaps
        BLOCK = sig > threshSIG;
        % remove the non-gaps that are less than 70
        BLOCK = bwareaopen(BLOCK,70);
        % close the pot holder chunks
        BLOCK = imclose(BLOCK,strel('disk',100));
        eBLOCK = BLOCK;
        %eBLOCK = imdilate(BLOCK,strel('disk',EXT));
        % make an image mask
        MASK = ~repmat(eBLOCK,[size(I,1) 1]);
        % get the bounding boxes for each mask
        R = regionprops(~MASK,'BoundingBox','PixelIdxList');
        for e = 1:numel(R)
            tmpMASK = zeros(size(MASK));
            tmpMASK(R(e).PixelIdxList) = 1;
            tmpMASK = imdilate(tmpMASK,strel('disk',EXT));
            tmpMASK(:,1:BKBOUND) = 0;
            tmpMASK(:,end-BKBOUND) = 0;
            tmpR = regionprops(logical(tmpMASK),'BoundingBox');
            R(e).BoundingBox = tmpR(1).BoundingBox;
        end
        fprintf(['starting: find the cone-tainer :' num2str(numel(R)) '\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find the cone-tainers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for each image block
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmpMASK_final = zeros(size(G));
        out = double(I)/255;
        HEIGHT = NaN*ones(1,numel(R));
        WIDTH = HEIGHT;
        dBIOMASS = WIDTH;
        plantHEIGHT = HEIGHT;
        % for each cone-tainer
        for e = 1:numel(R)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % crop a vertical strip for each container
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf(['starting: crop the strip\n']);
            % crop the strip
            tmpD = imcrop(I,R(e).BoundingBox);
            % next boundingbox
            nR(e).BoundingBox = R(e).BoundingBox;
            % trim the top
            tmpD(1:topTRIM:end,:,:) = [];
            fprintf(['ending: crop the strip\n']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % crop a vertical strip for each container
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % find the top of the cone-tainer - rough pass
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            BOTTOM_CUT_OFF = 200;
            fprintf(['starting: find the top of the cone-tainer\n']);
            MASK = getMASK_ver0(tmpD);
            MASK = imclose(MASK,strel('disk',5));
            fidx = find(MASK(end,:));
            MASK(end,fidx(1):fidx(end)) = 1;
            MASK = imfill(MASK,'holes');
            E = edge(MASK);
            sig = sum(E,2);
            sig = imfilter(sig,fspecial('average',[5 1]),'replicate');
            sig((end-BOTTOM_CUT_OFF):end) = 0;
            [J,nidx] = max(sig);
            nR(e).BoundingBox(4) = (nidx-OFFSET)-nR(e).BoundingBox(2);
            fprintf(['ending: find the top of the cone-tainer\n']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % find the top of the cone-tainer - rough pass
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % find plant thresholded on diameter of stem - fine pass
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf(['starting: get the plant mask\n']);
            cnt = 1;
            diameterThreshold = 40;
            STEM_SNIP = 20;
            outLN = 4;
            MAXcnt = 150; 
            SNIPBasePoint = 100;
            tmpD = imcrop(I,nR(e).BoundingBox);
            MASK = getMASK_ver0(tmpD);
            MASK = connectPlant(MASK);
            %MASK = bwlarge(MASK);
            [skeleton] = getPlantSkelton(MASK);
            [basePoint] = getStemBasePoint(MASK,SNIPBasePoint,skeleton);
            [diameter] = measureStemDiameter(MASK,STEM_SNIP,outLN,basePoint);
            fprintf(['clipping: the bottom plant mask\n']);
            while diameter > diameterThreshold && cnt < MAXcnt
                fprintf(['.']);
                [skeleton] = getPlantSkelton(MASK);
                [basePoint] = getStemBasePoint(MASK,SNIPBasePoint,skeleton);
                [diameter] = measureStemDiameter(MASK,STEM_SNIP,outLN,basePoint);
                if diameter > diameterThreshold
                    MASK(end,:) = [];
                    MASK = bwlarge(MASK);
                    nR(e).BoundingBox(4) = nR(e).BoundingBox(4) -1;
                end
                cnt = cnt + 1;
            end
            fprintf(['\n']);
            fprintf(['starting: get the plant mask\n']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % find plant thresholded on diameter of stem - fine pass
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % re-crop the plant,make mask,skeleton,basePoint,diameter and connectPlant
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tmpD = imcrop(I,nR(e).BoundingBox);
            MASK = getMASK_ver0(tmpD);
            MASK = connectPlant(MASK);
            MASK = bwlarge(MASK);
            [skeleton] = getPlantSkelton(MASK);
            [basePoint] = getStemBasePoint(MASK,SNIP,skeleton);
            [diameter] = measureStemDiameter(MASK,STEM_SNIP,outLN,basePoint);
            fprintf(['ending: get the plant mask\n']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % re-crop the plant,make mask,skeleton,basePoint,diameter and connectPlant
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % spool out the stem diameter
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            pNameDIA = [oPath nm '{PlantNumber_' num2str(e) '}{Phenotype_StemDiameter}.csv'];
            pushList{end+1} = pNameDIA;
            csvwrite(pNameDIA,diameter);
             
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % spool out the stem diameter
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % make phenotype measurements
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [plantHEIGHT(e),HEIGHT(e), WIDTH(e),dBIOMASS(e),hMassDistribution,vMassDistribution,CenterOfMass,StdDis] = measureHeightWidthDigitalBioMass(MASK);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % make phenotype measurements
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % spool out the stem diameter
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tmpName = [oPath nm '{PlantNumber_' num2str(e) '}{Phenotype_Std1}.csv'];
            csvwrite(tmpName,StdDis(1));
            pushList{end+1} = tmpName;
            tmpName = [oPath nm '{PlantNumber_' num2str(e) '}{Phenotype_Std2}.csv'];
            csvwrite(tmpName,StdDis(2));
            pushList{end+1} = tmpName;
            tmpName = [oPath nm '{PlantNumber_' num2str(e) '}{Phenotype_CenterOfMass1}.csv'];
            csvwrite(tmpName,CenterOfMass(1));
            pushList{end+1} = tmpName;
            tmpName = [oPath nm '{PlantNumber_' num2str(e) '}{Phenotype_CenterOfMass2}.csv'];
            csvwrite(tmpName,CenterOfMass(2));
            pushList{end+1} = tmpName;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % spool out the stem diameter
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % find the top of the container
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data.isPlant(e) = 0;
            if sum(MASK(:))/prod(size(MASK)) < thresP & dBIOMASS(e) > 100 & plantHEIGHT(e) > 20
                data.isPlant(e) = 1;
             
                skeleton = bwmorph(skeleton,'skel');
                for r = 1:15
                    ep = bwmorph(skeleton,'endpoints');
                    ep(end,:) = 0;
                    skeleton = skeleton - ep;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % trace the skeleton
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [path,pathcost,midx,DP] = traceSeedlingSkeleton(skeleton,basePoint);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % trace the skeleton
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                distT = bwdist(~MASK);
                
                pathStack = [];
                for r = 1:numel(path)
                    pathStack = [pathStack;path{r}];
                end
                UQ = unique(pathStack);
                pathCount = [];
                for u = 1:numel(UQ)
                    pathCount(u) = sum(pathStack==UQ(u));
                end
                fidx = find(pathCount == 1);
                
                
                UQu = UQ(fidx);
                pathStackU = [];
                for r = 1:numel(path)
                    tmpPath = path{r};
                    for u = 1:numel(UQu)
                        fidx = find(tmpPath == UQu(u));
                        pathStackU = [pathStackU;tmpPath(fidx)];
                    end
                end
                
                
                fidx = find(pathCount ~= 1);
                UQn = UQ(fidx);
                pathStackN = [];
                for r = 1:numel(path)
                    tmpPath = path{r};
                    for u = 1:numel(UQn)
                        fidx = find(tmpPath == UQn(u));
                        pathStackN = [pathStackN;tmpPath(fidx)];
                    end
                end
                
                
                leafOnlySkeleton = zeros(size(MASK));
                leafOnlySkeletonN = zeros(size(MASK));
                sidxU = sub2ind(size(MASK),DP(1,pathStackU)',DP(2,pathStackU)');
                sidxN = sub2ind(size(MASK),DP(1,pathStackN)',DP(2,pathStackN)');
                
                leafOnlySkeleton(sidxU) = 1;
                leafOnlySkeletonN(sidxN) = 1;
                leafOnlySkeletonN = imdilate(leafOnlySkeletonN,strel('disk',3));
                leafOnlySkeleton = logical((leafOnlySkeleton - leafOnlySkeletonN)==1);
                
                
                
                
                leafOnlySkeleton = bwareaopen(leafOnlySkeleton,50);
                leafOnlySkeleton = bwmorph(leafOnlySkeleton,'skel');
                for r = 1:20
                    ep = bwmorph(leafOnlySkeleton,'endpoints');
                    leafOnlySkeleton = leafOnlySkeleton - ep;
                end
                lR = regionprops(logical(leafOnlySkeleton),'PixelIdxList');
                midxMAX = [];
                for r = 1:numel(lR)
                    tmp = distT(lR(r).PixelIdxList);
                    [~,midxMAX(r)] = max(tmp);
                    midxMAX(r) = lR(r).PixelIdxList(midxMAX(r));
                end
                
                
                
                
                sidx = sub2ind(size(MASK),DP(1,pathStack)',DP(2,pathStack)');
                
                %sidx = find(skeleton);
                leafWidth = distT(sidx);
                rm = leafWidth < diameter/2;
                sidx(rm) = [];
                [l1 l2] = ind2sub(size(MASK),sidx);
                leafWidth(rm) = [];
                pathStack(rm) = [];
                [maxWLeafWidth,MAXIDX] = max(leafWidth);
                
                
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % make the mask overlay
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf(['starting: mask overlay\n']);
                % make mask overlay
                [x1 x2] = find(MASK);
                x1 = x1 + floor(nR(e).BoundingBox(2));
                x2 = x2 + floor(nR(e).BoundingBox(1));
                tmpMASK = zeros(size(tmpMASK_final));
                for i = 1:numel(x1)
                    tmpMASK(round(x1(i)),round(x2(i))) = 1;
                end
                tmpMASK_final = tmpMASK + tmpMASK_final;
                fprintf(['ending: mask overlay\n']);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % make the mask overlay
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % store the traced paths and other data
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf(['starting: storing the phenotypic data\n']);
                data.K{e} = [];
                for i = 1:numel(pathcost)
                    out = cwtK_imfilter(DP(:,path{i})',{5});
                    data.K{e} = [data.K{e};out.K];
                    dL = diff(DP(:,path{i}),1,2);
                    dL = sum(sum(dL.*dL,1).^.5);
                    data.pathLength{e}(i) = dL;
                    data.PATHS{e}{i} = [DP(2,path{i})+nR(e).BoundingBox(1);DP(1,path{i})+nR(e).BoundingBox(2)];
                end
                data.longestPathLength{e} = data.pathLength{e}(midx);
                data.longestPath(e) = midx;
                fprintf(['ending: storing the phenotypic data\n']);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % store the traced paths and other data
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % display the results for each plant
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf(['starting: display single plant\n']);
                yOFFSET = HEIGHT(e);
                yOFFSET = 1;
                xOFFSET = nR(e).BoundingBox(1);
                dE = bwboundaries(MASK);
                for r = 1:numel(dE)
                    tmp(r) = size(dE{r},1);
                end
                [J,sidx] = max(tmp);
               
                E = edge(MASK(yOFFSET:end,:));
                %E = imdilate(E,strel('disk',3,0));
                %out = flattenMaskOverlay(tmpD(yOFFSET:end,:,:), logical(MASK(yOFFSET:end,:)),.15,'g');
                %out = flattenMaskOverlay(out, logical(E),.55,'b');
                image(tmpD); hold on;
                for sidx = 1:numel(dE)
                    plot(dE{sidx}(:,2),dE{sidx}(:,1),'b')
                end
                
                axis off
                axis equal
                hold on
                Xbar = nR(e).BoundingBox(1):nR(e).BoundingBox(1)+nR(e).BoundingBox(3) - xOFFSET;
                Ybar = (nR(e).BoundingBox(2)+nR(e).BoundingBox(4))*ones(size(Xbar)) - yOFFSET;

                %{
                % something is not right here
                plot(Xbar,Ybar,'r')
                Ybar = (HEIGHT(e))*ones(size(Xbar));
                plot(Xbar,Ybar,'r')
                %}

                for i = 1:numel(data.PATHS{e})
                    plot(data.PATHS{e}{i}(1,:) - xOFFSET,data.PATHS{e}{i}(2,:) - yOFFSET,'k','LineWidth',1);
                end
                plot(data.PATHS{e}{data.longestPath(e)}(1,:)-xOFFSET,data.PATHS{e}{data.longestPath(e)}(2,:)-yOFFSET,'r','LineWidth',1);
                title(['plant ' num2str(e)]);
                set(gca,'Position',[0 0 1 1]);
                plot(CenterOfMass(2),CenterOfMass(1),'c*');
                plot(linspace((CenterOfMass(2)-WIDTH(e)/2),(CenterOfMass(2)+WIDTH(e)/2),2),[CenterOfMass(1) CenterOfMass(1)],'c');
                plot([CenterOfMass(2) CenterOfMass(2)],linspace((CenterOfMass(1)-plantHEIGHT(e)/2),(CenterOfMass(1)+plantHEIGHT(e)/2),2),'c');
                TH = linspace(-pi,pi,200);
                X = StdDis(2)*cos(TH) + CenterOfMass(2);
                Y = StdDis(1)*sin(TH) + CenterOfMass(1);
                plot(X,Y,'c')
                plot(basePoint(2),basePoint(1)-20,'c.');
                plot(linspace((basePoint(2)-diameter/2),(basePoint(2)+diameter/2),2),[basePoint(1)-20 basePoint(1)-20],'c');
                
                plot(l2(MAXIDX),l1(MAXIDX),'g*');
                
                %plot(DP(2,pathStack),DP(1,pathStack),'k.');
                
                for r = 1:numel(midxMAX)
                    [subIDX1 subIDX2] = ind2sub(size(MASK),midxMAX(r));
                    plot(subIDX2,subIDX1,'g*')
                end
                
                plot(DP(2,pathStack(MAXIDX)),DP(1,pathStack(MAXIDX)),'g*');
                tmpImageName = [oPath nm '{PlantNumber_' num2str(e) '}{Phenotype_Image}.tif'];
                
                saveas(gca,tmpImageName);
                pushList{end+1} = tmpImageName;
                close all
                fprintf(['ending: display single plant\n']);
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % display the results - for each plant
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
        
    catch ME
        getReport(ME)
        close all
    end
    
    
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % display the results rendered onto the raw image
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close all
        fprintf(['starting with image and results display \n']);
        out = flattenMaskOverlay(I, logical(tmpMASK_final),.15,'g');
        h = image(out);
        axis off
        hold on
        imwrite(tmpMASK_final,[oPath nm '{PlantNumber_All}{Phenotype_ImageMask}.tif']);
        for e = 1:numel(HEIGHT)
            if ~(plantHEIGHT(e)==0)
                % make the data for the bars for height
                Xbar = nR(e).BoundingBox(1):nR(e).BoundingBox(1)+nR(e).BoundingBox(3);
                Ybar = (nR(e).BoundingBox(2)+nR(e).BoundingBox(4))*ones(size(Xbar));
                % plot bars
                plot(Xbar,Ybar,'r')
                Ybar = (HEIGHT(e))*ones(size(Xbar));
                plot(Xbar,Ybar,'r')
                
                if data.isPlant(e)
                    % plot paths for tracing
                    for i = 1:numel(data.PATHS{e})
                        plot(data.PATHS{e}{i}(1,:),data.PATHS{e}{i}(2,:),'r','LineWidth',2);
                    end
                    % plot longest path
                    plot(data.PATHS{e}{data.longestPath(e)}(1,:),data.PATHS{e}{data.longestPath(e)}(2,:),'k','LineWidth',2);
                end
                fprintf(['ending with image and results display \n']);
            end
        end
        axis equal;axis off;drawnow;set(gca,'Position',[0 0 1 1]);
        tmpImageName = [oPath nm '{PlantNumber_All}{Phenotype_Image}.tif'];
        saveas(gca,tmpImageName);
        pushList{end+1} = tmpImageName;
        close all
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % display the results
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pName1 = [oPath nm '{PlantNumber_All}{Phenotype_DigitalBioMass}.csv'];
        pName2 = [oPath nm '{PlantNumber_All}{Phenotype_PlantHeight}.csv'];
        pName3 = [oPath nm '{PlantNumber_All}{Phenotype_LongestPath}.csv'];
        pNameWidth = [oPath nm '{PlantNumber_All}{Phenotype_Width}.csv'];
        csvwrite(pName1,dBIOMASS);
        pushList{end+1} = pName1;
        csvwrite(pName2,plantHEIGHT);
        pushList{end+1} = pName2;
        csvwrite(pName3,data.longestPathLength);
        pushList{end+1} = pName3;
        csvwrite(pNameWidth,WIDTH);
        pushList{end+1} = pNameWidth;
        for e = 1:numel(HEIGHT)
            if (plantHEIGHT(e)==0)
                data.K{e} = NaN;
            end
            pName4 = [oPath nm '{PlantNumber_' num2str(e) '}{Phenotype_Curvature}.csv'];
            csvwrite(pName4,data.K{e});
            pushList{end+1} = pName4;
        end
        pName5 = [oPath nm '{PlantNumber_All}{Phenotype_imageFileName}.txt'];
        fileID = fopen(pName5,'w');
        nbytes = fprintf(fileID,'%s\n',imageFile);
        pushList{end+1} = pName5;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % push to iRODS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pushToiRods(rPath,pushList);
    catch ME
        getReport(ME)
        close all
    end
    close all
end


%{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compile
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    compile_directory = '/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/Maize/maizeSeedling/tmpSubmitFiles/';
    CMD = ['mcc -d ' compile_directory ' -a im2single.m -m -v -R -singleCompThread singleSeedlingImage.m'];
    eval(CMD);


    oPath = '/mnt/spaldingdata/nate/mirror_images/maizeData/hirsc213/return/seedlingData/output/';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/growthchamber12.7.15/day11_b73_10-11.nef';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC4.10.16/plot100.nef';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC3.5.16/plot19c.tiff';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC3.8.16/plot19c.tiff';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC3.11.16/plot19c.tiff';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC4.5.16/plot100.nef';
    %fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC4.6.16/plot100.nef';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC4.2.16/plot100.nef';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC4.8.16/plot100.nef';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC3.28.16/plot100.nef';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC2.24.16/plot9c.tiff';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/growthchamber12.7.15/day12_ph207_1-2-3comp.tif';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC3.17.16/plot54c.tiff'; % bad qr
    fileName ='/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC3.12.16/plot40c.tiff';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC2.25.16/plot9c.tiff';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC2.27.16/plot10c.tiff';
    filename = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC2.27.16/plot11c.tiff';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC3.2.16/plot19c.tiff';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC2.26.16/plot14c.tiff';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC2.22.16/plot11c.tiff';
    %fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC2.29.16/plot23c.tiff';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC2.26.16/plot7c.tiff';
    %fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC5.6.16/plot224.nef';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC2.26.16/plot10c.tiff';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC4.10.16/plot105.nef';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC4.1.16/plot108.nef';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC4.1.16/plot110.nef';
    %fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC4.1.16/plot112.nef';
    %fileName = '/iplant/home/hirsc213/maizeData/seedlingData/one_month_test_mac/GC2.26.16/plot1.nef';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/13-Jul-2016/{Plot_388}{Experiment_16}{Planted_7-1-16}{SeedSource_DI1911-9}{SeedYear_2015}{Genotype_PH207}{Treatment_Control}{PictureDay_13}.nef'
    

    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/01-Aug-2016/{Plot_401}{Experiment_17}{Planted_7-15-16}{SeedSource_26 DI1926 x Mo17}{SeedYear_2015}{Genotype_B73 x Mo17}{Treatment_Control}{PictureDay_17}.nef';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/05-Aug-2016/{Plot_401}{Experiment_17}{Planted_7-15-16}{SeedSource_26 DI1926 x Mo17}{SeedYear_2015}{Genotype_B73 x Mo17}{Treatment_Control}{PictureDay_21}.nef';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/04-Aug-2016/{Plot_392}{Experiment_17}{Planted_7-15-16}{SeedSource_DI1909-11}{SeedYear_2015}{Genotype_B73}{Treatment_Control}{PictureDay_20}.nef';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/22-Jul-2016/{Plot_392}{Experiment_17}{Planted_7-15-16}{SeedSource_DI1909-11}{SeedYear_2015}{Genotype_B73}{Treatment_Control}{PictureDay_7}.nef';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/30-Jul-2016/{Plot_401}{Experiment_17}{Planted_7-15-16}{SeedSource_26 DI1926 x Mo17}{SeedYear_2015}{Genotype_B73 x Mo17}{Treatment_Control}{PictureDay_15}.nef';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/31-Jul-2016/{Plot_393}{Experiment_17}{Planted_7-15-16}{SeedSource_DI1909-11}{SeedYear_2015}{Genotype_B73}{Treatment_Control}{PictureDay_16}.nef';
    singleSeedlingImage(fileName,100,5,400,100,4,20,'/mnt/spaldingdata/nate/seedlingData/','/iplant/home/hirsc213/maizeData');


    FilePath = '/mnt/spaldingdata/nate/mirror_images/maizeData/hirsc213/return/seedlingData/output/';
    txtFileList = {};
    FileExt = {'txt'};
    txtFileList = gdig(FilePath,txtFileList,FileExt,1);
    [imageFile] = getImageFile(txtFileList,'10','17');
    [imageFile] = getImageFile(txtFileList,'10','19');
    [imageFile] = getImageFile(txtFileList,'10','9'); % not present
    [imageFile] = getImageFile(txtFileList,'106','9');
    [imageFile] = getImageFile(txtFileList,'100','22');
    [imageFile] = getImageFile(txtFileList,'224','21');
    [imageFile] = getImageFile(txtFileList,'104','19'); % for 18 - cant read QR code
    [imageFile] = getImageFile(txtFileList,'104','12'); % for 11 - cant read
    [imageFile] = getImageFile(txtFileList,'105','21');
    [imageFile] = getImageFile(txtFileList,'108','12');
    [imageFile] = getImageFile(txtFileList,'110','12'); % wrong read pDay11
    [imageFile] = getImageFile(txtFileList,'112','12');
    [imageFile] = getImageFile(txtFileList,'1','11');
    [imageFile] = getImageFile(txtFileList,'100','14');
    [imageFile] = getImageFile(txtFileList,'100','7');
    [imageFile] = getImageFile(txtFileList,'102','18');
    singleSeedlingImage(imageFile,100,5,100,100,4,oPath);
    parfor e = 1:22
        [imageFile] = getImageFile(txtFileList,'100',num2str(e));
        if ~isempty(imageFile)
            singleSeedlingImage(imageFile,100,5,100,100,4,oPath);
        end
    end
    
%}