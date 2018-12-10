function [phenotypes,measureBlocks,JSONdoc] = measurePhenotypes(MASK,SKELETON,I,LEN)
    SNIP = 10;
    STEM_SNIP = 10;
    outLN = 4;

    % setup measureBlocks
    measureBlocks(1).history = [0 1 2];
    measureBlocks(1).statBlock = [zeros(8,8,8)];

    measureBlocks(2).history = [0 2 4];
    measureBlocks(2).statBlock = [zeros(8,8,8)];

    measureBlocks(3).history = [0 4 8];
    measureBlocks(3).statBlock = [zeros(8,8,8)];

    
    JSONdoc = [];
    phenotypes = [];

    for e = 1:numel(MASK)
        try
            % get the distance transform
            DIST = double(bwdist(~logical(MASK{e})));
            % get the basepoint
            [basePoint] = getStemBasePoint(MASK{e},SNIP,SKELETON{e});
            % get the stem diameter
            [diameter] = measureStemDiameter(MASK{e},STEM_SNIP,outLN,basePoint);
            % trace the skeleton
            [path,pathcost,midx,DP] = traceSeedlingSkeleton(SKELETON{e},basePoint);
            for p = 1:numel(path)
                colorFunc = [];
                tmpy = DP(:,path{p})';
                for k = 1:size(I{e},3)
                    colorFunc(:,k) = interp2(I{e}(:,:,k),tmpy(:,2),tmpy(:,1));
                end
                % measure kurvature
                out = cwtK_imfilter(DP(:,path{p})',{5});
                DIS = interp2(DIST,DP(2,path{p})',DP(1,path{p})');
                K{p} = [DP(:,path{p})' out.K DIS colorFunc];
            end
            
            [mC pidx] = max(pathcost);
            maxPath = DP(:,path{pidx});
            dP = diff(maxPath,1,2);
            L = sum(sum(dP.*dP,1).^.5);
           
            
            
            pathValues = interp2(double(DIST),maxPath(2,:)',maxPath(1,:)');
            
            
            % extract color information
            HSV = rgb2hsv(I{e});
            RGB = I{e};
            LAB = rgb2lab(RGB);
            tmpMask = MASK{e};
            fidx = find(find(tmpMask));

            rgbVec = [];
            hsvVec = [];
            labVec = [];
            for k = 1:3
                tmpi = RGB(:,:,k);
                rgbVec = [rgbVec;tmpi(fidx)'];

                tmpi = HSV(:,:,k);
                hsvVec = [hsvVec;tmpi(fidx)'];

                tmpi = LAB(:,:,k);
                labVec = [labVec;tmpi(fidx)'];

            end

            histoRGB = [];
            histoHSV = [];
            histoLAB = [];
            for k = 1:3
                histoRGB(k,:) = hist(rgbVec(k,:),linspace(0,1,255));
                histoHSV(k,:) = hist(hsvVec(k,:),linspace(0,1,255));
                histoLAB(k,:) = hist(labVec(k,:),linspace(0,1,255));
            end

            RE = regionprops(tmpMask,'Perimeter','Solidity','MajorAxisLength','MinorAxisLength','EquivDiameter','Eccentricity');
            

            % measure some classic phenotypes
            [plantHEIGHT,HEIGHT, WIDTH,dBIOMASS,vMassDistribution,hMassDistribution,CenterOfMass,StdDis] = measureHeightWidthDigitalBioMass(MASK{e});
           

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tmpDoc = [];
            % register the object as a plant
            tmpDoc = generatePhenotypeNode(tmpDoc,1,{'isPlant'},'isPlant');
            % save the base point
            tmpDoc = generatePhenotypeNode(tmpDoc,basePoint,{'plantBasePoint'},'basePoint');
            % save the stem diameter
            tmpDoc = generatePhenotypeNode(tmpDoc,diameter,{'stemDiameter'},'stemDiameter');
            % save the plant height
            tmpDoc = generatePhenotypeNode(tmpDoc,plantHEIGHT,{'plantHeight'},'plantHeight');
            % save the plant width
            tmpDoc = generatePhenotypeNode(tmpDoc,WIDTH,{'plantWidth'},'plantWidth');
            % save the digital biomass
            tmpDoc = generatePhenotypeNode(tmpDoc,dBIOMASS,{'digitalBiomass'},'digitalBiomass');
            % save the horizontal mass distribution
            newHD = interp1(1:numel(hMassDistribution),hMassDistribution,linspace(1,numel(hMassDistribution),1000));
            tmpDoc = generatePhenotypeNode(tmpDoc,newHD,{'horizontalMassDistribution'},'horizontalMassDistribution');
            % save the vertical mass distribution
            newVD = interp1(1:numel(vMassDistribution),vMassDistribution,linspace(1,numel(vMassDistribution),1000));
            tmpDoc = generatePhenotypeNode(tmpDoc,newVD,{'verticalMassDistribution'},'verticalMassDistribution');
            % save the center of mass
            tmpDoc = generatePhenotypeNode(tmpDoc,CenterOfMass,{'centerOfMass'},'centerOfMass');
            % save the standard dev of mass
            tmpDoc = generatePhenotypeNode(tmpDoc,StdDis,{'standardDeviation'},'standardDeviation');
            % save the standard dev of mass
            tmpDoc = generatePhenotypeNode(tmpDoc,L,{'longestPathLength'},'longestPathLength');
            % save the color data - rgb
            tmpDoc = generatePhenotypeNode(tmpDoc,histoRGB',{'RGBcolorValue','channel'},'RGBcolorValue');
            % save the color data - histoHSV
            tmpDoc = generatePhenotypeNode(tmpDoc,histoHSV',{'HSVcolorValue','channel'},'HSVcolorValue');
            % save the color data - histoHSV
            tmpDoc = generatePhenotypeNode(tmpDoc,LEN',{'ResolutionConversion'},'ResolutionConversion');
            % render the classic binary morphology
            flds = fields(RE);
            for f = 1:numel(flds)
                tmpDoc = generatePhenotypeNode(tmpDoc,RE.(flds{f}),{flds{f}},flds{f});
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [JSONdoc] = appendStruct(JSONdoc,tmpDoc);
            %JSONdoc(e) = tmpDoc;


        

            % store phenotypes
            phenotypes(e).isPlant = 1;
            phenotypes(e).basePoint = basePoint;
            phenotypes(e).stemDiameter = diameter;
            phenotypes(e).plantHEIGHT = plantHEIGHT;
            phenotypes(e).plantWIDTH = WIDTH;
            phenotypes(e).digitalBIOMASS = dBIOMASS;
            phenotypes(e).hMassDistribution = hMassDistribution;
            phenotypes(e).vMassDistribution = vMassDistribution;
            phenotypes(e).CenterOfMass = CenterOfMass;
            phenotypes(e).StdDis = StdDis;
            phenotypes(e).longestPathLength = L;
            phenotypes(e).longestPath = maxPath;
            phenotypes(e).distanceTransformAlongPath = pathValues;
            phenotypes(e).histoRGB = histoRGB;
            phenotypes(e).histoHSV = histoHSV;
            phenotypes(e).histoLAB = histoLAB;
            flds = fields(RE);
            for f = 1:numel(flds)
                phenotypes(e).(flds{f}) = RE.(flds{f});
            end

            % stack all paths
            for p = 1:numel(path)
                %phenotypes(e).paths(p).d = DP(:,path{p});
                phenotypes(e).paths(p).d = K{p};
            end

            failedList(e) = 0;
        catch ME
            tmpDoc = [];
            % register the object as a plant
            tmpDoc = generatePhenotypeNode(tmpDoc,0,{'isPlant'},'isPlant');

            [JSONdoc] = appendStruct(JSONdoc,tmpDoc);
            %JSONdoc(e) = tmpDoc;


            failedList(e) = 1;
            getReport(ME)
%{
            junk.isPlant = 0;
            junk
            [phenotypes] = appendStruct(phenotypes,junk);
%}
        end
    end
    
    for e = 1:numel(failedList)
        phenotypes(e).isPlant = ~logical(failedList(e));
    end
    
     
end


%{

            %{
            disp = 0;
            for p = 1:numel(path)
                mpath = convertPath(DP(:,path{p}));
                [measureBlocks prob(p)] = measurePathStats(fliplr(mpath),measureBlocks,toM);
                
                if disp
                    imshow(SKELETON{e},[])
                    hold on
                    plot(DP(2,path{p}),DP(1,path{p}),'r')
                    title(num2str(prob(p)))
                    drawnow
                    waitforbuttonpress
                end
                
            end
            %}
            
            
            
            %{
            [DP,T] = makeSkeletonTraceValues(SKELETON{e});
            
            
            [value,sidx,path] = pathSnap(T,sourceIDX,targetsIDX);
            
            
            
            
            [endPoints(:,1),endPoints(:,2)] = find(bwmorph(SKELETON{e},'endpoints'));
            
            
            
            
            
            [branchPoints(:,1),branchPoints(:,2)] = find(bwmorph(SKELETON{e},'branchpoints'));
            %}
            

%{
            %bwmorph(SKELETON{e},'spur')
            
            % get the frequency use of points
            pH = zeros(1,size(DP,2));
            for p = 1:numel(path)
                pH(path{p}) = pH(path{p}) + 1;
            end
            
            
            
            % find the points that are used most frequently
            fidx = find(pH > 2);
            
            
            % calculate the stalklocation
            stalkLocation = mean(DP(:,fidx),2);
            stalkHeight = max(DP(:,fidx),[],2);
            
            
            % calculate the number of branches, branch angle, and leaf
            % curvature
            [branchPoints(:,1),branchPoints(:,2)] = find(bwmorph(SKELETON{e},'branchpoints'));
            
            stalk = [stalkLocation(2)*ones(size(MASK{e},1),1),(1:size(MASK{e},1))'];
            
            
            
            for p = 1:size(branchPoints,1)
                tmpD = bsxfun(@minus,stalk,flip(branchPoints(p,:),2));
                tmpD = (sum(tmpD.*tmpD,2)).^.5;
                [dis,midx] = min(tmpD);
                if dis < 10
                    
                end
                
            end
            
            
            
            
            
            
            
            imshow(SKELETON{e},[]);
            hold on
            plot(stalkLocation(2)*ones(1,size(MASK{e},1)),1:size(MASK{e},1),'r');
            %}

%}